# Copyright (c) 2018 Robert Mieth and Yury Dvorkin
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# +++++
# tools.jl
#
# Some auxillary useful functions
#


# Empirical Residual Distribution
# Input:    residuals ... matrix (n x t)
#           bus ... scalar : number of bus
# Output:   F, domain ...  F: vector of values F(s), domain: vector s
function get_F_for_bus(residuals, bus)
    res_at_bus = ϵ_hat[bus,:]
    min_res = minimum(res_at_bus)
    max_res = maximum(res_at_bus)
    steps = 100
    stepsize = abs(max_res - min_res)/100
    F = []
    domain = collect(min_res:stepsize:max_res)
    for s in domain
        F_s = count(x -> x <= s, res_at_bus)
        push!(F, F_s)
    end
    F = 1/T .* F
    return F, domain
end

# Quantile estimation
function get_quantile_for_bus(η, residuals, bus)
    res_at_bus = ϵ_hat[bus,:]
    res_sorted = sort(res_at_bus)
    l = length(res_sorted)
    quantile = res_sorted[ceil(Int, η*l)]
    return quantile
end

#Build Distance Matrix
function get_distance_matrix(n_buses, lines)
    G = Graph(n_buses)
    for l in values(lines)
        add_edge!(G, l.to_node, l.from_node)
    end
    dijkstra_shortest_paths(G, 1).dists
    dist_matrix = []
    for v in vertices(G)
        append!(dist_matrix, gdistances(G, v))
    end
    dist_matrix = reshape(dist_matrix, n_buses, n_buses)
    return dist_matrix
end

# Multivariaate distribution can not handle zero variance
# Zero has to be added after generating the errors
function get_noise(mvdist, zero_indices)
    noise = rand(mvdist)
    final_length = length(noise) + length(zero_indices)
    v = zeros(final_length)
    skip = 0
    for i in 1:final_length
        if i in zero_indices
            v[i] = 0
            skip += 1
        else
            v[i] = noise[i-skip]
        end
    end
    return v
end

function create_random_init_data(t_init, no_load_buses, mvdist; lambda_max=100)
    x_hist = []
    λ_hist = []
    for t in 1:t_init
        #Some random prices
        λ_t = [rand()*lambda_max for i in 1:n_buses]
        x_t = x_det(λ_t) + get_noise(mvdist, no_load_buses)
        if t==1
            x_hist = x_t
            λ_hist = λ_t
        else
            x_hist = hcat(x_hist, x_t)
            λ_hist = hcat(λ_hist, λ_t)
        end
    end
    return λ_hist, x_hist
end
