using CSV
using DataFrames
using Graphs
using XLSX


# get the DiscreteGraphAlgorithms package - REQUIRED TO USE THESE FUNCTIONS
using DiscreteGraphAlgorithms





##################################################################
#    NOTE: START FUNCTIONS, BASED ON GRAPHS.JL IMPLEMENTATION    #
##################################################################

function _accumulate_basic_subgroup!(
    betweenness::Vector{Float64}, 
    state::Graphs.DijkstraState, 
    g::AbstractGraph,
    si::Integer,
    vs::Vector{Int64},
)
    n_v = length(state.parents) # this is the ttl number of vertices
    δ = zeros(n_v)
    σ = state.pathcounts
    P = state.predecessors

    # make sure the source index has no parents.
    P[si] = []
    # we need to order the source vertices by decreasing distance for this to work.
    S = reverse(state.closest_vertices) # Replaced sortperm with this
    for w in S
        
        !(w in vs) && continue
            
        coeff = (1.0 + δ[w]) / σ[w]
        for v in P[w]
            if v > 0
                δ[v] += (σ[v] * coeff)
            end
        end
        if w != si
            betweenness[w] += δ[w]
        end
    end
    return nothing
end



function betweenness_centrality_subgroup(
    g::AbstractGraph,
    vs = vertices(g),
    distmx::AbstractMatrix=weights(g);
    normalize=true,
    #endpoints=false,
)
    vs = collect(vs)
    n_v = nv(g)
    k = length(vs)
    isdir = is_directed(g)

    betweenness = zeros(n_v)
    for s in vs
        if degree(g, s) > 0  # this might be 1?
            state = dijkstra_shortest_paths(g, s, distmx; allpaths=true, trackvertices=true)
            #if endpoints
            #    _accumulate_endpoints!(betweenness, state, g, s)
            #else
            _accumulate_basic_subgroup!(betweenness, state, g, s, vs)
            #end
        end
    end

    _rescale_subgroup!(betweenness, n_v, normalize, isdir, k, vs, )

    return betweenness
end



function build_centralities_subgroup(
    df::DataFrame,
    graph_wrapper::GraphWrapper;
    field_key::Symbol = :key,
    field_subtype::Symbol = :sub_type_1,
)
    
    graph = graph_wrapper.graph
    df_base = filter(x -> !ismissing(x[field_subtype]), df)
    all_subtypes = collect(sort(unique(df_base[:, field_subtype])))

    # try converting
    try 
        all_subtypes = String.(all_subtypes)
    catch e
        @info "Error converting type to String in build_centralities_subgroup(): trying to conver to Symbol, then String..."
        all_subtypes = String.(Symbol.(all_subtypes))
        @info "Success."
    end
    
    # set matrix
    new_mat = zeros(Float64, graph_wrapper.dims[1], length(all_subtypes))
    
    # iterate over available subtypes
    for (j, subtype) in enumerate(all_subtypes)
        
        verts_group = get_verts_by_type(
            df_base,
            graph_wrapper,
            subtype;
            field = field_subtype,
        )

        vec = betweenness_centrality_subgroup(
            graph,
            verts_group;
            normalize = false,
        )
        
        new_mat[:, j] = vec # ordered by graph
    end
    
    # build matrix out
    new_mat = DataFrame(new_mat, Symbol.(all_subtypes))
    new_mat[:, field_key] = graph_wrapper.vertex_names
    
    df_base = leftjoin(
        df_base,
        new_mat,
        on = [field_key]
    )
    
    return df_base
end



function _rescale_subgroup!(
    betweenness::Vector{Float64}, 
    n::Integer, 
    normalize::Bool, 
    directed::Bool, 
    k::Integer,
    vs::Vector{Int64},
)
    if normalize
        if n <= 2
            do_scale = false
        else
            do_scale = true
            scale_in = 1.0 / ((k - 1) * (k - 2))
            scale_out = 1.0 / (k * (k - 1))
        end
    else
        if !directed
            do_scale = true
            scale = 1.0 / 2.0
        else
            do_scale = false
        end
    end
    if do_scale
        #if k > 0
        #    scale = scale * n / k
        #end
        
        v_scale_in = Int64.(in.(1:n_v, (vs, )))
        v_scale_out = (v_scale_in .+ 1) .% 2
        scale = scale_in * v_scale_in + scale_out * v_scale_out
        
        betweenness .*= scale
    end
    
    return nothing
end



function get_verts_by_type(
    df::DataFrame,
    graph_wrapper::GraphWrapper,
    subtype::String;
    field::Symbol = :sub_type_1,
)
    # get 
    group = filter(
        x -> (x[field] == subtype),
        df
    )[:, :key]

    verts = findall(
        x -> (x in group),
        graph_wrapper.vertex_names,
    )
    
    verts_complement = setdiff(vertices(graph_wrapper.graph), verts);
    
    return verts_complement
    
end






"""


# Constructs

```
write_subgroup_centralities(
    fp_in::String,
    fp_out::String,
    sheet_name_vertices::String,
    sheet_name_graph::String;
    field_key::Symbol = :key,
    field_subtype::Symbol = :sub_type_1,
)
```

##  Function Arguments

- `fp_in`: input path to Excel file storing data
- `fp_out`: output path to CSV storing summary dataframe
- `sheet_name_graph`: sheet name in `fp_in` with graph sparse adjacency
- `sheet_name_vertices`: sheet name in `fp_in` with vertex information


##  Keyword Arguments

- `field_key`: field in vertex table (`sheet_name_vertices`) that stores the key 
    for the vertex (name)
- `field_source`: field in graph table (`sheet_name_graph`) that stores the source
    along edges in the sparse adjacency matrix
- `field_subtype`: field in vertex table (`sheet_name_vertices`) that stores the 
    vertex type. Types are used to group vertices and calculate subgroup 
    centralities
- `field_target`: field in graph table (`sheet_name_graph`) that stores the target
    along edges in the sparse adjacency matrix
- `return_df`: specify as true to return the output dataframe
"""
function write_subgroup_centralities!(
    fp_in::String,
    fp_out::String,
    sheet_name_graph::String, # "2010-2019A"
    sheet_name_vertices::String; # "Nodes"
    field_key::Symbol = :key,
    field_source::Symbol = :source,
    field_subtype::Symbol = :sub_type_1,
    field_target::Symbol = :target,
    return_df::Bool = false,
)
    
    # get graph and vertices
    df_graph_cur = XLSX.readtable(fp_in, sheet_name_graph) |> DataFrame;
    df_vertices = XLSX.readtable(fp_in, sheet_name_vertices) |> DataFrame;
    
    # build graph
    gw_full = DiscreteGraphAlgorithms.df_to_graph_wrapper(
        df_graph_cur, 
        field_source, 
        field_target,
    );
    
    # iterate to build subgroup centralities
    df_out = build_centralities_subgroup(
        df_vertices,
        gw_full;
        field_key = field_key,
        field_subtype = field_subtype,
    )
    
    # write to output CSV
    CSV.write(fp_out, df_out)
    
    
    out = return_df ? df_out : nothing
    
    return out
end

