# Copyright (c) 2018 Robert Mieth and Yury Dvorkin
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# +++++
# testcase.jl
#
# Default testcase specification
# Invoke file using include()

function return_case_data()
    # Name of the case
    case_id = "testcase"
    exp_id = "testcase"

    # Specify data files
    datadir  = "data/feeder_data/basecase_lv_noneg"
    price_file = "data/price_data/rand_max200_min30_n10000.csv"

    # total number of timesteps
    t_total = 50
    # number of inital timesteps (>=2)
    t_init = 2

    # Model settings
    robust_cc = true # If true run using der dist. robust OPF, else derterministic
    enable_voltage_constraints = true # if false, voltage constraints are never enforced
    enable_generation_constraints = true # if false only det. generation constraints are enforced (if not generation constraints are enforced the problem is unbounded)
    enable_flow_constraints = true # if false, flow constraints are never enforced
    compare_to_detopf = false # only DR signals are derived from the robust model, generation dispatch is again taken from a det. OPF run
    run_power_flow_test = true # check result feasibility after runs

    # Voltage at root bus
    v_root = 1
    # Costumer Tariff
    tariff = 25
    # Price for demand Response
    dr_price = 150
    # Inital dr price assumption
    dr_price_assumption = 150
    # Voltage Security margin
    η_v = 0.1
    # Generation Security margin
    η_g = 0.1
    # Demand standard deviation relative to load
    relative_std = 0.1
    # Partizipation Factor
    α = zeros(15)
    α[1] = 1
    # Correlation settings
    max_correlation = 0
    # Factor for higher load
    load_fact = 1

    # Create vector of parameters for each bus
    # Change here for individual prices at each bus
    β1_set = ones(15)./dr_price
    β0_set = zeros(15)

    β1_init = ones(15)./dr_price_assumption
    β0_init = zeros(15)

    return case_id, exp_id, datadir, price_file, t_total, t_init, robust_cc, enable_voltage_constraints, enable_generation_constraints, enable_flow_constraints, compare_to_detopf, run_power_flow_test, v_root, tariff, η_v, η_g, relative_std, α, max_correlation, load_fact, β1_set, β0_set, β1_init, β0_init
end

return_case_data()
