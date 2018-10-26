<<<<<<< HEAD
# invoke using include()

function return_case_data()
    # Name of the case
    case_id = "only_voltage"
    exp_id = "fullrun1_d2"

    # Specify data files
    datadir  = "data/feeder_data/basecase_lv_noneg"
    price_file = "data/price_data/rand_max200_min30_n10000.csv"

    # total number of timesteps
    t_total = 500
    # number of inital timesteps (>=2)
    t_init = 2

    # Model settings
    robust_cc = true
    enable_voltage_constraints = true
    enable_generation_constraints = false
    enable_flow_constraints = false
    compare_to_detopf = true
    run_power_flow_test = true

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
    relative_std = 0.2
    # Partizipation Factor
    α = zeros(15)
    α[1] = 1
    # Correlation settings
    max_correlation = 0
    # Factor for higher load
    load_fact = 2

    β1_set = ones(15)./dr_price
    β0_set = zeros(15)

    β1_init = ones(15)./dr_price_assumption
    β0_init = zeros(15)

    return case_id, exp_id, datadir, price_file, t_total, t_init, robust_cc, enable_voltage_constraints, enable_generation_constraints, enable_flow_constraints, compare_to_detopf, run_power_flow_test, v_root, tariff, η_v, η_g, relative_std, α, max_correlation, load_fact, β1_set, β0_set, β1_init, β0_init
end

return_case_data()
