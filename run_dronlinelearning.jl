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
using Distributions # For error generating Distribution
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

# Result storages for regression
x_hist = Dict()
λ_hist = Dict()
β1_hat_hist = Dict()
β0_hat_hist = Dict()
# Result storages for OPF results
results_hist = []
results_oracle_hist = []
results_varorc_hist = []
results_psorc_hist = []
status_hist = DataFrame(t=[], orac=[], varorc=[], psorc=[], learning=[])
# Result storages for PF results
outcome_hist = []
outcome_oracle_hist = []
outcome_varorc_hist = []
outcome_psorc_hist = []

Σ_oracle = Σ_true
μ_oracle = μ_true
Ω_oracle = Ω_true

@assert length(β0_init) == n_buses
@assert length(β1_init) == n_buses

println(">>>>> Starting loop for $(t_total) timestep$(t_total>1?"s":"")")
for t in 1:t_total
    println(">>>>> Running for t=$t (of $t_total)")

    # Regression phase
    β0_hat_t = Dict("learning" => zeros(n_buses), "varorc" => zeros(n_buses), "psorc" => β0, "oracle" => β0)
    β1_hat_t = Dict("learning" => zeros(n_buses), "varorc" => zeros(n_buses), "psorc" => β0, "oracle" => β0)
    # Only 'learning' and 'varorc' need to estimate the prices
    for info in ["learning", "varorc"]
        if t <= t_init || info ∉ keys(λ_hist)
            # Use assumption for inital timesteps
            β0_hat_t[info] = β0_init
            β1_hat_t[info] = β1_init
        else
            # Perform Regression with the historical data
            for i in 1:n_buses
                # Matrix implementation
                # fisher = [sum(λ_hist[info][i,k]^2 for k in 1:T)  sum(λ_hist[info][i,:]); sum(λ_hist[info][i,:]) T]
                # obs = [sum(λ_hist[info][i,:] .* x_hist[info][i,:]) ; sum(x_hist[info][i,:])]
                # est = (fisher^-1)*obs
                # β1_hat_t[info][i] = est[1]
                # β0_hat_t[info][i] = est[2]

                # Explicit implementation (generally faster)
                λ_mean = mean(λ_hist[info][i,:])
                x_mean = mean(x_hist[info][i,:])
                β1_hat_t[info][i] = sum((λ_hist[info][i,:] .- λ_mean) .* (x_hist[info][i,:] .- x_mean))/sum((λ_hist[info][i,:] .- λ_mean).^2)
                β0_hat_t[info][i] = x_mean - β1_hat_t[info][i]*λ_mean

                # if regression was not viable in t use estimate from last iteration
                isnan(β1_hat_t[info][i]) && (β1_hat_t[info][i] = β1_hat_hist[info][i,t-1])
                isnan(β0_hat_t[info][i]) && (β1_hat_t[info][i] = β0_hat_hist[info][i,t-1])
            end
        end
        # Truncate
        β1_min = 1/10000
        β1_hat_t[info] = max.(β1_hat_t[info], β1_min)
        β0_hat_t[info] = max.(β1_hat_t[info], 0)

        # Store
        if t == 1
            β1_hat_hist[info] = β1_hat_t[info]
            β0_hat_hist[info] = β0_hat_t[info]
        else
            temp = pop!(β1_hat_hist, info)
            β1_hat_hist[info] = hcat(temp, β1_hat_t[info])
            temp = pop!(β0_hat_hist, info)
            β0_hat_hist[info] = hcat(temp, β0_hat_t[info])
        end
    end

    # Store oracle values
    for info in ["oracle", "psorc"]
        if t == 1
            β1_hat_hist[info] = β1_hat_t[info]
            β0_hat_hist[info] = β0_hat_t[info]
        else
            temp = pop!(β1_hat_hist, info)
            β1_hat_hist[info] = hcat(temp, β1_hat_t[info])
            temp = pop!(β0_hat_hist, info)
            β0_hat_hist[info] = hcat(temp, β0_hat_t[info])
        end
    end

    # Estimate Error variance
    μ_hat = Dict("learning" => zeros(n_buses), "oracle" => zeros(n_buses), "varorc" => zeros(n_buses), "psorc" => zeros(n_buses))
    Σ_hat = Dict("learning" => zeros(n_buses, n_buses), "oracle" => zeros(n_buses, n_buses), "varorc" => zeros(n_buses, n_buses), "psorc" => zeros(n_buses, n_buses))
    Ω_hat = Dict("learning" => zeros(n_buses+1, n_buses+1), "oracle" => zeros(n_buses+1, n_buses+1), "varorc" => zeros(n_buses+1, n_buses+1), "psorc" => zeros(n_buses+1, n_buses+1))
    # Only 'learning' and 'psorc' need to estimate the covariance matrix
    for info in ["learning", "psorc"]
        if t <= 2 || info ∉ keys(λ_hist)
            # No residuals in first round
            μ_hat[info] = zeros(n_buses)
            Σ_hat[info] = zeros(n_buses, n_buses)
        else
            # x_hats resulting from parameter estimates
            x_hat = []
            available_samples = size(λ_hist[info], 2)
            for s in 1:available_samples
                x_hat_t = β1_hat_t[info] .* λ_hist[info][:,s] .+ β0_hat_t[info]
                if s == 1
                    x_hat = x_hat_t
                else
                    x_hat = hcat(x_hat, x_hat_t)
                end
            end

            # Resulting residuals
            ϵ_hat = x_hist[info] - x_hat
            # ϵ_hat = map(x -> abs(x)>num_thershold ? x : 0., ϵ_hat)   # Get rid of numerical noise
            μ_hat[info] = [mean(ϵ_hat[i,:]) for i in 1:n_buses]
            # μ_hat = map(x -> abs(x)>num_thershold ? x : 0., μ_hat) # Get rid of numerical noise
            # var_e_hat = [var(ϵ_hat[i,:]) for i in 1:n_buses]
            # var_e_hat = map(x -> abs(x)>num_thershold ? x : 0., var_e_hat)
            Σ_hat[info] = cov(ϵ_hat') # Emp Covarince
        end
        # Second order moment matrix of residuals
        Ω_hat[info] = [Σ_hat[info]+μ_hat[info]*μ_hat[info]' μ_hat[info]; μ_hat[info]' 1] # Second order moment Matrix
    end

    # Get wholesale price
    wholesale = ws_prices[t]

    # Run the OPF models
    result_learning, status_learning, solvetime_learning = run_demand_response_opf(feeder, β1_hat_t["learning"], β0_hat_t["learning"],
                            wholesale, μ_hat["learning"], Σ_hat["learning"], Ω_hat["learning"];
                            α=α, x_in=[], model_type="x_opf", robust_cc=false, enable_flow_constraints=true,
                            enable_voltage_constraints=true, enable_generation_constraints=true)

    result_oracle, status_oracle, solvetime_oracle = run_demand_response_opf(feeder, β1, β0,
                            wholesale, μ_oracle, Σ_oracle, Ω_oracle;
                            α=α, x_in=[], model_type="x_opf", robust_cc=false, enable_flow_constraints=true,
                            enable_voltage_constraints=true, enable_generation_constraints=true)


    display(result_oracle)

    # Get and safe the price results
    λ_t_learning = max.(result_learning[:lambda],0)
    if t == 1
        λ_hist["learning"] = λ_t_learning
    else
        temp = pop!(λ_hist, "learning")
        λ_hist["learning"] = hcat(temp, λ_t_learning)
    end
    # λ_t_oracle = results_oracle[:lambda]
    # λ_t_oracle = map(x -> x < 0 ? 0 : x, λ_t_oracle)

    # Observe an reaction (so that error is the same for all)
    ϵ_t = get_noise(ϵ_mvdist, no_load_buses)
    x_obs_learning = x_det(λ_t_learning) .+ ϵ_t
    # x_obs_oracle = x_det(λ_t_oracle) .+ ϵ_t
    # x_obs_varorc = x_det(λ_t_varorc) .+ ϵ_t
    # x_obs_psorc = x_det(λ_t_psorc) .+ ϵ_t

    if t == 1
        x_hist["learning"] = x_obs_learning
    else
        temp = pop!(x_hist, "learning")
        x_hist["learning"] = hcat(temp, x_obs_learning)
    end
# repeat
end


d = Dict("a"=> [1,2], "b" => 1)
pop!(d, "b")
d["b"] = "hallo"
