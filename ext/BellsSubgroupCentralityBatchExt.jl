"""
Calculate Bells Subgroup Centralities in batch using an input excel. 
"""
module BellsSubgroupCentralityBatchExt


using CSV
using DataFrames
using Graphs
using XLSX
# get the DiscreteGraphAlgorithms package - REQUIRED TO USE THESE FUNCTIONS
# this is why this is an optional extension
using DiscreteGraphAlgorithms

#path = dirname(@__FILE__)
#include(joinpath(path, "bells_centrality.jl"))


export write_subgroup_centralities!


##################################################################
#    NOTE: START FUNCTIONS, BASED ON GRAPHS.JL IMPLEMENTATION    #
##################################################################

#=
EXPLORED APPROACH FROM GRAPHS.jl TO MODIFY
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
=#

#=
EXPLORED APPROACH FROM GRAPHS.jl TO MODIFY

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
=#



"""
Function underlying the 
"""
function build_centralities_subgroup(
    df::DataFrame,
    graph_wrapper::GraphWrapper;
    field_key::Symbol = :key,
    field_subtype::Symbol = :sub_type_1,
)
    
    graph = graph_wrapper.graph

    df_base = copy(df)
    filter!(x -> !ismissing(x[field_subtype]), df_base)
   
    # try converting types
    try 
        df_base[!, field_subtype] = String.(df_base[:, field_subtype])

    catch e
        @info "Error converting type to String in build_centralities_subgroup(): trying to convert to Symbol, then String..."

        df_base[!, field_subtype] = String.(
            Symbol.(
                df_base[:, field_subtype]
            )
        )

        @info "Success."
    end

    all_subtypes = collect(sort(unique(df_base[:, field_subtype])))
    
    # set matrix
    n_col = length(all_subtypes)*3 + 1
    new_mat = zeros(Float64, graph_wrapper.dims[1], n_col)
    new_mat_fields = Vector{String}(["" for x in 1:n_col])
    new_mat_fields[1] = "betweenness_centrality"

    # calculate betweenness centrality for each vertex
    vec_betweenness = betweenness_centrality(
        graph;
        normalize = false,
    )
    new_mat[:, 1] = vec_betweenness

    # iterate over available subtypes
    for (j, subtype) in enumerate(all_subtypes)
        
        # get matrix indices and fields
        inds = 3*j - 2 .+ collect(1:3)
        fields = build_subtype_names(subtype)

        verts_group, verts_complement = get_verts_by_type(
            df_base,
            graph_wrapper,
            subtype;
            field = field_subtype,
        )
        
        # get local centrality (B_s) 
        vec_local = betweenness_centrality_bells(
            graph;
            vertices_source = verts_group,
            vertices_target = verts_group,
        )

        # get global centrality (B_s^C) 
        vec_global = betweenness_centrality_bells(
            graph;
            vertices_source = verts_complement,
            vertices_target = verts_complement,
        )

        # get boundary centrality (Bb_s) 
        vec_boundary = betweenness_centrality_bells(
            graph;
            vertices_source = verts_group,
            vertices_target = verts_complement,
        )

        vec_boundary += betweenness_centrality_bells(
            graph;
            vertices_source = verts_complement,
            vertices_target = verts_group,
        )
        
        # local, global, and boundary
        new_mat[:, inds[1]] = vec_local # ordered by graph
        new_mat[:, inds[2]] = vec_global # ordered by graph
        new_mat[:, inds[3]] = vec_boundary # ordered by graph
        
        new_mat_fields[inds[1]:inds[3]] = fields
    end
    
    # build matrix out
    new_mat = DataFrame(new_mat, Symbol.(new_mat_fields))
    new_mat[:, field_key] = graph_wrapper.vertex_names
    
    df_base = leftjoin(
        df_base,
        new_mat,
        on = [field_key]
    )
    
    return df_base
end



function build_subtype_names(
    subtype::String,
)
    out = [
        "local_bs",
        "global_bsc",
        "boundary_bbs"
    ]

    out = ["$(subtype)_$(x)" for x in out]

    return out
end


#=
FROM GRAPHS.jl -- EXPLORED MODIFYING IT

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
        
        v_scale_in = Int64.(in.(1:n, (vs, )))
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
    out = (verts, verts_complement)

    return out
end
=#





"""
Read an input Excel file, which splits vertices and edges into different sheets,
    and calculate subgroup centralities across different groups (which can be 
    specified using field_subtype, which defaults to `field_subtype_1`)


# Constructs

```
write_subgroup_centralities!(
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
- `field_source`: field in graph table (`sheet_name_graph`) that stores the 
    source along edges in the sparse adjacency matrix
- `field_subtype`: field in vertex table (`sheet_name_vertices`) that stores the 
    vertex type. Types are used to group vertices and calculate subgroup 
    centralities
- `field_target`: field in graph table (`sheet_name_graph`) that stores the 
    target along edges in the sparse adjacency matrix
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


end
