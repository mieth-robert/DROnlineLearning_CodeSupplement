# OnlineLearningDRSupplement
Respostitory containing supplementary data and code for "Online Learning for Network Constrained Demand Response Pricing in Distribution Systems" by Mieth and Dvorkin. [DOI](pending)

*Installation instructions:*

The optimization model was implemented by using [JuMP](https://github.com/JuliaOpt/JuMP.jl) and auxillary packages in the [Julia](http://julialang.org/downloads/) programming language.
Additionally, we used Mosek 8.1.0.63 in our numerical experiments. [Mosek](https://www.mosek.com) is a commercial solver which must be installed and licensed. The solver was choosen for its specifc features for semi-definite programming. For more information on solvers, see the JuMP documentation.

The experiments requre Julia 0.6 and the following Julia packages:
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

* Running the code:*

