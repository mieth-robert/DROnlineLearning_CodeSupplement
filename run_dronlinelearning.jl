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
include("src/test_power_flow.jl")

if length(ARGS) == 0 println(">>>>> No case argument provided. Proceeding with default testcase") end

#1 Load case settings
casefile = length(ARGS) > 0 ? ARGS[1] : "cases/testcase.jl"

case_id, exp_id, datadir, price_file, t_total, t_init, robust_cc, enable_voltage_constraints,
        enable_generation_constraints, enable_flow_constraints, compare_to_detopf, run_power_flow_test, v_root, tariff, η_v,
        η_g, relative_std, α, max_correlation, β1_set, β0_set, β1_init, β0_init = include(casefile)
println(">>>>> Succesfully loaded $(case_id)")

#2 Load Data from file
# Feeder Data
feeder = load_feeder(datadir)
# Price Data
ws_prices, timestamps = read_price_data(price_file)

#3 Prepare everything for simulation
t_total == 0 && (t_total = length(ws_prices))

# Feeder Data
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
result_hist = Dict()
status_hist = DataFrame(t=[], learning=[], oracle=[], varorc=[], psorc=[])
solvetime_hist = DataFrame(t=[], learning=[], oracle=[], varorc=[], psorc=[])
# Result storages for PF results
outcome_hist = Dict()
var_hat_hist = []

# Create files for results storage
results_path = "results/$(exp_id == "" ? Dates.format(Dates.now(), "mm-dd_HHMM") : exp_id)/$case_id"
mkpath(results_path)

Σ_oracle = Σ_true
μ_oracle = μ_true
Ω_oracle = Ω_true

@assert length(β0_init) == n_buses
@assert length(β1_init) == n_buses

println(">>>>> Starting loop for $(t_total) timestep$(t_total>1?"s":"")")
tic()
for t in 1:t_total
    println(">>>>> Running for t=$t (of $t_total)")

    # Regression phase
    β0_hat_t = Dict("learning" => zeros(n_buses), "varorc" => zeros(n_buses), "psorc" => β0, "oracle" => β0)
    β1_hat_t = Dict("learning" => zeros(n_buses), "varorc" => zeros(n_buses), "psorc" => β1, "oracle" => β1)
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
                isnan(β0_hat_t[info][i]) && (β0_hat_t[info][i] = β0_hat_hist[info][i,t-1])
            end
        end
        # Truncate
        β1_min = 1/10000
        β1_hat_t[info] = max.(β1_hat_t[info], β1_min)
        # β0_hat_t[info] = max.(β0_hat_t[info], 0)

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
    set_root_feeder_cost(feeder, wholesale)

    # Run the OPF models
    result = Dict(); status = Dict(); solvetime = Dict()
    println(">>>>> Starting Learning Model")
    result["learning"], status["learning"], solvetime["learning"] = run_demand_response_opf(feeder, β1_hat_t["learning"], β0_hat_t["learning"],
                            μ_hat["learning"], Σ_hat["learning"], Ω_hat["learning"];
                            α=α, x_in=[], model_type="x_opf", robust_cc=robust_cc, enable_flow_constraints=enable_flow_constraints,
                            enable_voltage_constraints=enable_voltage_constraints, enable_generation_constraints=enable_generation_constraints)

    println(">>>>> Starting Oracle Model")
    result["oracle"], status["oracle"], solvetime["oracle"] = run_demand_response_opf(feeder, β1, β0,
                            μ_oracle, Σ_oracle, Ω_oracle;
                            α=α, x_in=[], model_type="x_opf", robust_cc=robust_cc, enable_flow_constraints=enable_flow_constraints,
                            enable_voltage_constraints=enable_voltage_constraints, enable_generation_constraints=enable_generation_constraints)

    println(">>>>> Starting Ps-oracle Model")
    result["psorc"], status["psorc"], solvetime["psorc"] = run_demand_response_opf(feeder, β1, β0,
                            μ_hat["psorc"], Σ_hat["psorc"], Ω_hat["psorc"];
                            α=α, x_in=[], model_type="x_opf", robust_cc=robust_cc, enable_flow_constraints=enable_flow_constraints,
                            enable_voltage_constraints=enable_voltage_constraints, enable_generation_constraints=enable_generation_constraints)

    println(">>>>> Starting Var-oracle Model")
    result["varorc"], status["varorc"], solvetime["varorc"] = run_demand_response_opf(feeder, β1_hat_t["varorc"], β0_hat_t["varorc"],
                            μ_oracle, Σ_oracle, Ω_oracle;
                            α=α, x_in=[], model_type="x_opf", robust_cc=robust_cc, enable_flow_constraints=enable_flow_constraints,
                            enable_voltage_constraints=enable_voltage_constraints, enable_generation_constraints=enable_generation_constraints)

    # Create vector of errors that is the same for all cases
    ϵ_t = get_noise(ϵ_mvdist, no_load_buses)
    for info in ["learning", "oracle", "psorc", "varorc"]
        info in keys(result) || continue
        
        # Determine system settings from det OPF
        if compare_to_detopf
            result_detopf, status_detopf, solvetime_detopf = run_demand_response_opf(feeder, β1, β0,
                                μ_oracle, Σ_oracle, Ω_oracle;
                                α=α, x_in=result[info][:x_opt], model_type="opf", robust_cc=false, enable_flow_constraints=true,
                                enable_voltage_constraints=true, enable_generation_constraints=true)
        end

         # Get and save the price results
        λ_t = max.(result[info][:lambda],0)
        if t == 1
            λ_hist[info] = λ_t
        else
            temp = pop!(λ_hist, info)
            λ_hist[info] = hcat(temp, λ_t)
        end
    
        # Observe an reaction disturbed by the error vector
        x_obs = x_det(λ_t) .+ ϵ_t

        # Save the reactions
        if t == 1
            x_hist[info] = x_obs
        else
            temp = pop!(x_hist, info)
            x_hist[info] = hcat(temp, x_obs)
        end

        # Run power flow using the true reactions
        outcome = run_power_flow_dr(feeder, (compare_to_detopf ? result_detopf : result[info]), x_obs)

        # Save to history
        result[info][:t] = fill(t, nrow(result[info]))
        outcome[:t] = fill(t, nrow(outcome))
        if t==1 
            result_hist[info] = result[info]
            outcome_hist[info] = outcome
        else
            result_hist[info] = vcat(result_hist[info], result[info])
            outcome_hist[info] = vcat(outcome_hist[info], outcome)
        end
        # Save to hdd
        # CSV.write("$(results_path)/results_$(info).csv", result[info], append=!(t==1))
        # CSV.write("$(results_path)/outcomes_$(info).csv", outcome, append=!(t==1))

    end
    push!(status_hist, [t, status["learning"], status["oracle"], status["varorc"], status["psorc"]])
    push!(solvetime_hist, [t, solvetime["learning"], solvetime["oracle"], solvetime["varorc"], solvetime["psorc"]])

# repeat
end
toc()
display(status_hist)

# save to hdd
println(">>>>> Saving results")
for info in ["learning", "oracle", "psorc", "varorc"]
    CSV.write("$(results_path)/results_$(info).csv", result_hist[info])
    CSV.write("$(results_path)/outcomes_$(info).csv", outcome_hist[info])
    beta_0_hist_df = DataFrame()
    beta_1_hist_df = DataFrame()
    for b in 1:n_buses
        beta_0_hist_df[Symbol("bus_$b")] = β0_hat_hist[info][b,:]
        beta_1_hist_df[Symbol("bus_$b")] = β1_hat_hist[info][b,:]
    end
    CSV.write("$(results_path)/beta_1_hat_$(info).csv", beta_1_hist_df)
    CSV.write("$(results_path)/beta_0_hat_$(info).csv", beta_0_hist_df)
end  

CSV.write("$(results_path)/status.csv", status_hist)
CSV.write("$(results_path)/solvetime.csv", solvetime_hist)
wholesale_df = DataFrame(wholesale=ws_prices[1:t_total])
CSV.write("$(results_path)/wholesale_prices.csv", wholesale_df)
