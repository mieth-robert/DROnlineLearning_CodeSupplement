# OnlineLearningDRSupplement
Repository containing supplementary data and code for "Online Learning for Network Constrained Demand Response Pricing in Distribution Systems" by Robert Mieth and Yury Dvorkin.

*Installation instructions:*

The optimization model was implemented by using [JuMP](https://github.com/JuliaOpt/JuMP.jl) and auxiliary packages in the [Julia](http://julialang.org/downloads/) programming language.
Additionally, we used Mosek 8.1.0.63 in our numerical experiments. [Mosek](https://www.mosek.com) is a commercial solver which must be installed and licensed. The solver was chosen for its specific features for semi-definite programming. For more information on solvers, see the JuMP documentation.

The experiments require Julia 0.6 and the following Julia packages:
- [JuMP](https://github.com/JuliaOpt/JuMP.jl)0.18.3
- [Mosek](https://github.com/JuliaOpt/Mosek.jl)0.8.4
- [Distributions](https://github.com/JuliaStats/Distributions.jl)0.15.0
- [CSV](https://github.com/JuliaData/CSV.jl)0.2.5
- [DataFrames](https://github.com/JuliaData/DataFrames.jl)0.11.7

You should force the use of particular versions of these Julia packages with 
```
julia> Pkg.pin("JuMP", v"0.18.3")
julia> Pkg.pin("Mosek", v"0.8.4")
julia> Pkg.pin("Distributions", v"0.15.0")
julia> Pkg.pin("CSV", v"0.2.5")
julia> Pkg.pin("DataFrames", v"0.11.7")
```

*Running the code:*


The code for the experiments is contained in the ``run_dronlinelearning.jl`` file as well as the auxiliary files in the ``src`` directory. 
The ``run_dronlinelearning.jl`` file is the main file that loads the data and performs the iterative algorithm based on the setting specified in a case file. The file  ``input.jl`` defines necessary data types and functions to load an prepare the network data, ``model_definitions.jl`` contains the JuMP formulations, ``test_power_flow.jl`` contains the methods for feasibility checks and ``tools.jl`` contains some auxiliary functions.

You can run a testcase provided in ``cases/testcase.jl`` by executing:
```
julia run_dronlinelearning.jl cases/testcase.jl
```
If no specific case argument is provided, the script will run the testcase per default. The casefile is also a julia file and is heavily commented. It allows to easily parametrize the experiment. Note, that a single timestep may take up to several seconds in computation. 

