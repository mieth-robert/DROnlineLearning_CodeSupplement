# Copyright (c) 2018 Robert Mieth and Yury Dvorkin
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# +++++
# run_dronlinelearning.jl
#
# This model refers to the paper: XXX
#
# +++++
# Devnotes: ---

# Invocing Julia Packages
using JuMP  # Optimization Package
using Mosek # Solver Interface, requires Mosek to be installed
using CSV, DataFrames   # For Handling Data
using Distributions
using LightGraphs # Utilizing some predefined algorithms

include("src/input.jl")
include("src/model_definitions.jl")
include("src/tools.jl")

if length(ARGS) == 0 println(">>>>> No case argument provided. Proceeding with default testcase") end

#1 Load case settings
casefile = length(ARGS) > 0 ? ARGS[1] : "cases/testcase.jl"

case_id, datadir, price_file, t_total, t_init, robust_cc, enable_voltage_constraints,
        enable_generation_constraints, enable_flow_constraints, run_power_flow_test, v_root, tariff, η_v,
        η_g, relative_std, α, max_correlation, β1_set, β0_set, β1_init, β0_init = include(casefile)
println(">>>>> Succesfully loaded $(case_id)")

#2 Load Data from file
# Feeder Data
feeder = load_feeder(datadir)
# Price Data
ws_prices, timestamps = read_price_data(price_file)

#3 Prepare everything for simulation
buses = feeder.buses
lines = feeder.lines
n_buses = feeder.n_buses

# Find nodes without load (hence no demand response)
loads = [b.d_P for b in buses]
non_dr_buses = []
no_load_buses = []
for (i,b) in enumerate(buses)
    abs(b.d_P) > 0 && continue
    push!(non_dr_buses, i)
    push!(no_load_buses, i)
    # Those two are probably the same but not necessarily; I'll leave it like this
    # if there is need for differnetiation later
end

# Check and set demand response behavior
@assert length(β1_set) == n_buses
@assert length(β0_set) == n_buses
β1 = β1_set
β0 = β0_set
std_vec = loads .* relative_std

# Define true cov matrix and noise generating distribution
dist_matrix = get_distance_matrix(n_buses, lines)
max_dist = maximum(dist_matrix)
corr_step = max_correlation/max_dist
correlation_factors = dist_matrix.*(max_correlation/max_dist)
Σ_true = zeros(n_buses, n_buses)
μ_true = zeros(n_buses)
r_matrix = zeros(n_buses, n_buses)
for i in 1:n_buses
    for j in 1:n_buses
        if i == j
            r = 1
        else
            r = (max_dist - dist_matrix[i,j] + 1)/max_dist * corr_step
        end
        covij = std_vec[i]*std_vec[j]*r
        Σ_true[i,j] = covij
        Σ_true[j,i] = covij
        r_matrix[i,j] = r
        r_matrix[j,i] = r
    end
end
Σ_reduced = Σ_true[setdiff(1:end, no_load_buses), setdiff(1:end, no_load_buses)]
ϵ_mvdist = MvNormal(zeros(n_buses-length(no_load_buses)), Σ_reduced)
# Define response generating function
x_det(λ) = (β1 .* λ) .+ β0
# Second order moment matrix
Ω_true = [Σ_true+μ_true*μ_true' μ_true; μ_true' 1]

# Result storages
x_hist = []
λ_hist = []
results_hist = []
results_oracle_hist = []
results_varorc_hist = []
results_psorc_hist = []
lines_hist = []
lines_oracle_hist = []
β1_hat_hist = []
β0_hat_hist = []
status_hist = DataFrame(t=[], orac=[], varorc=[], psorc=[], learning=[])
outcome_hist = []
outcome_oracle_hist = []
outcome_varorc_hist = []
outcome_psorc_hist = []

Σ_oracle = Σ_true
μ_oracle = μ_true
Ω_oracle = Ω_true
wholesale = 80


println(">>>>> Starting loop for $(t_total) timestep$(t_total>1?"s":"")")
for t in 1:t_total

    x = zeros(n_buses)
    result, status, solvetime = run_demand_response_opf(feeder, β1, β0, wholesale, μ_oracle, Σ_oracle, Ω_true;
                            α=α, x_in=x, model_type="opf", robust_cc=false, enable_flow_constraints=true,
                            enable_voltage_constraints=true, enable_generation_constraints=true)
    display(result)

end