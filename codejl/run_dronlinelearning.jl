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

cd("codejl")

# Invocing Julia Packages
using JuMP  # Optimization Package
using Mosek # Solver Interface, requires Mosek to be installed
using CSV, DataFrames   # For Handling Data

include("input.jl")

#1 Load case settings
    # TODO
datadir  = "../data/basecase_lv_noneg"

#2 Load Data from file
# Feeder Data
feeder = load_feeder(datadir)

# Price Data
