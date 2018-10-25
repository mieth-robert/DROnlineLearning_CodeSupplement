# invoke using include()

function return_case_data()
    # Name of the case
    case_id = "testcase"
    exp_id = "nydata"

    # Specify data files
    datadir  = "data/feeder_data/basecase_lv_noneg"
    price_file = "data/price_data/nyiso_nyc_2018-10-08_to_2018-10-10.csv"

    # total number of timesteps
    # set to 0 to use all availabl time steps from 'price_file'
    t_total = 0
    # number of inital timesteps (>=2)
    t_init = 2

    # Model settings
    robust_cc = true
    enable_voltage_constraints = true
    enable_generation_constraints = true
    enable_flow_constraints = true

    run_power_flow_test = true

    # Voltage at root bus
    v_root = 1
    # Costumer Tariff
    tariff = 25
    # Price for demand Response
    dr_price = 45
    # Inital dr price assumption
    dr_price_assumption = 45
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

    β1_set = ones(15)./dr_price
    β0_set = zeros(15)

    β1_init = ones(15)./dr_price_assumption
    β0_init = zeros(15)

    return case_id, exp_id, datadir, price_file, t_total, t_init, robust_cc, enable_voltage_constraints, enable_generation_constraints, enable_flow_constraints, run_power_flow_test, v_root, tariff, η_v, η_g, relative_std, α, max_correlation, β1_set, β0_set, β1_init, β0_init
end

return_case_data()
