
# invoke using include()

function return_case_data()
    # Name of the case
    case_id = "testcase"

    # Specify data files
    datadir  = "data/feeder_data/basecase_lv_noneg"
    price_file = "data/price_data/rand_max200_min30_n10000.csv"

    t_total = 10

    # Model settings

    # Increase load by this factor
    load_fact = 1

    # Voltage Security margin
    η_v = 0.1
    # Generation Security margin
    η_g = 0.1


    # Power System Settings
    v_root = 1

    return case_id, datadir, price_file, v_root
end

return_case_data()
