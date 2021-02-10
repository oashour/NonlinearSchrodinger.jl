# Visualization

## Recipes
`NonlinearSchrodinger.jl` provides recipes for `Plots.jl` to make plotting easy. You can simply call any plotting function on a `::Sim` or `::Calc` object, as shown in the [Examples](@ref) page.

## Visualizing Recursion

We can visualize the DT recursion as follows:

```@example1
using LightGraphs, MetaGraphs, Plots, GraphRecipes, LaTeXStrings
default(size=(1000,1000))

col = []

function calc_rs(n, p, g, N, pal, size, col)
    println("Calculating for ($n,$p) with order $N")
    if !haskey(g[:name], (n, p)) # if there is not already a node with name (n, p)
        # add a new vertex and set it's name to (n, p)
        # the most recent added vertex always has the index nv(g)
        add_vertex!(g) 
        set_prop!(g, nv(g), :name, (n, p))
        push!(size, N)
        push!(col, pal[N])
     end
     if n != 1
         calc_rs(n-1, 1, g, N-1, pal, size, col)
         calc_rs(n-1, p+1, g, N-1, pal, size, col)
         # at this point, the vertices with the names (n, p), (n-1, 1) and (n-1, p+1) already exist
         # so we look them up by their name
         v1 = g[(n, p), :name]
         v2 = g[(n-1, 1), :name]
         v3 = g[(n-1, p+1), :name]
               
         # add edges (n, p) -> (n-1, 1)  and (n, p) -> (n-1, p+1)
         add_edge!(g, v1, v2)
         add_edge!(g, v1, v3)
    end
    return nothing
end

N = 4
g = MetaDiGraph() # empty graph
set_indexing_prop!(g, :name) # allows one to look up nodes by the attribute :name
size = []
pal = palette(:heat, 7)
calc_rs(N, 1, g, N, pal, size, col)
```

```@example1
graphplot(g, names = [latexstring(get_prop(g, v, :name)) for v in vertices(g)], arrow=:arrow, 
              nodeshape=:circle, curvature_scalar=0.0, nodeweights=size, markercolor=col, 
              fontsize=14, markersize=0.04, linewidth=2, method=:tree, curves=false)
savefig("recurse_4.svg") #hide
```
![](recurse_4.svg)
