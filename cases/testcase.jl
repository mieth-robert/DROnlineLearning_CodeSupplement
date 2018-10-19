
# invoke using include()

function return_case_data()
    case_id = "testcase"

    datadir  = "data/feeder_data/basecase_lv_noneg"
    price_file = "data/price_data/rand_max200_min30_n10000.csv"

    return case_id, datadir, price_file
end

return_case_data()
