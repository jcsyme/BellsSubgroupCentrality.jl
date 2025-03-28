_DICT_VALID_SUBGROUP_TYPES = Dict(
    "boundary" => :b,
    "b" => :b,
    "global" => :g,
    "g" => :g,
    "local" => :l,
    "l" => :l,
    "overall" => :o,
    "o" => :o,
)



"""
Calculate the component betweenness centrality for paths starting at i and heading
    to all other vertices
"""
function b_i!(
    graph::AbstractGraph, # graph_wrapper_bell.graph
    vertex::Int64,
    vec_component::Vector{Float64};
    vertices_target::Union{Vector{Int64}, Nothing} = nothing,
)
    
    # initialize components
    
    # calculate shortest paths
    state = dijkstra_shortest_paths(
        graph, 
        vertex;
        allpaths = true,
    )
    
    vertices_target = isa(vertices_target, Nothing) ? collect(vertices(graph)) : vertices_target
    
    # iterate over predecessors 
    for j in vertices_target
        
        preds = state.predecessors[j]
        
        # skip if there is no path
        (length(preds) == 0) && continue
        
        vec_paths = Vector{Vector{Int64}}()
        vec_parents = Vector{Int64}()
        
        # get all paths
        paths_all = get_paths(
            state.predecessors, 
            preds,
            [j],
            vec_parents,
            vec_paths,
        )
        
        dict_counts = get_path_counts_by_vertex(vec_paths, )
        
        for (k, v) in dict_counts
            vec_component[k] += v/state.pathcounts[j]
        end
    end
    
    return nothing
end



"""
Recursively get paths 
# Constructs

get_paths(
    states::Vector, 
    node::Vector,
    parent::Vector{Int64},
    parents::Vector,
    paths::Vector{Int64},
)


##  Function Arguments

states: vector of states 

node: 

##  Keyword Arguments

"""
@inline function get_paths(
    states::Vector{Vector{Int64}}, 
    node::Vector{Int64},
    parent::Vector{Int64},
    parents::Vector{Int64},
    paths::Vector{Vector{Int64}},
)
    # here is where you want to check the node's name against your list of
    # 'nodes to shout out about'
    
    (length(node) == 0) && (return 1)
    
    
    @inbounds for (i, child) in enumerate(node)
        
        #push!(parents, parent...)
        parents = parent
        
        #print("i = $(i)\nchild = $(child)\nnext = $(states[child])\n")
        val = get_paths(
                states, 
                states[child], 
                [parent; [child]],
                parents,
                paths,
        )
            

        (val == 1) && (push!(paths, parents[2:end]); empty!(parents))

    end
    
    return 0
end



"""
Using a vector of vectors of ordered paths, calculate the number of times
    each vertex is on a path (keys in output dictionary)
"""
function get_path_counts_by_vertex(
    vec_paths::Vector,
)
    dict_out = Dict()
    
    # iter
    for (i, vec) in enumerate(vec_paths)
        for el in vec
            (el in keys(dict_out)) ? (dict_out[el] += 1) : (dict_out[el] = 1)
        end
    end
    
    return dict_out
end
    


"""
Calculate any one of the following subgroup centralities:

    * Overall (:overall or :o):
        Standard betweenness centrality
    * Local (:local or :l): 
        Betweenness centralities calculated only within the in-group
    * Global (:global or :g):
        Betweenness centralities calculated only within the out-group
    * Boundary subgroup centrality (:boundary or :b): 
        Betweenness centralities calculated between the in- and out-groups

See:

Bell, JR. Subgroup centrality measures. Network Science. 2014;2(2):277-297. 
doi:10.1017/nws.2014.15


# Construct

```
get_subgroup_centrality(
    graph::AbstractGraph,
    type::Symbol,
    vertices_in_group::Vector{Int64},
)
```

##  Function Arguments

- `graph`: AbstractGraph 
- `type`: one of the following types to calculate (not case sensitive)
    * :boundary or :b
    * :global or :g
    * :local or :l
    * :overall or :o
- `vertices_in_group`: vector of vertices that are considered in-group

##  Keyword Arguments

- `vertices_out_group`: optional passage of outgroup vertices. Otherwise, 
    calculated as vertex set complement of in-groups
"""
function get_subgroup_centrality(
    graph::AbstractGraph,
    type::Symbol,
    vertices_in_group::Vector{Int64};
    vertices_out_group::Union{Vector{Int64}, Nothing} = nothing,
)
    # get 
    isa(vertices_out_group, Nothing) && (vertices_out_group = setdiff(vertices(graph), vertices_in_group))
    type = get_subgroup_type(type)

    # case where the centrality keeps the same outputs
    if type in [:g, :l, :o]
        
        # switch on case (o -> nothing, g => out_group, l => in_group)
        verts = nothing
        if type in [:g, :l]
            verts = (type == :g) ? vertices_out_group : vertices_in_group
        end
        
        # call the function and return
        out = subgroup_centrality_bells(
            graph;
            vertices_source = verts,
            vertices_target = verts,
        )

        return out
    end
    

    ##  HERE'S THE CASE FOR BOUNDARY

    out = subgroup_centrality_bells(
        graph;
        vertices_source = vertices_in_group,
        vertices_target = vertices_out_group,
    )

    out .+= subgroup_centrality_bells(
        graph;
        vertices_source = vertices_out_group,
        vertices_target = vertices_in_group,
    )

    return out

end



"""
Check that the subgroup is specified correctly.
"""
function get_subgroup_type(
    type::Symbol,
)
    # check the type
    type = lowercase(String(type))
    out = get(_DICT_VALID_SUBGROUP_TYPES, type, nothing)

    if isa(out, Nothing) 
        msg = join(keys(_DICT_VALID_SUBGROUP_TYPES), ", ")
        error("Invalid type $(type) specified. Must be one of the following (case insensitive): $(msg).")
    end

    return out
end



"""
Explicit implementation of subgroup betweeness centrality that can be used to 
    implment Bell's centrality measures. Can be used to calculate the following 
    subgroup centralities:

    * Overall (0): 
        `vertices_source` and `vertices_target` are `nothing` (same as standard
        betweenness centrality)
    * Local (L): 
        `vertices_source` and `vertices_target` are each set to the in-group
    * Global (G):
        `vertices_source` and `vertices_target` are each set to the out-group 
        (vertices outside of the group)
    * Boundary subgroup centrality (B): 
        Sum of these two cases:
            `vertices_source = in_group`, `vertices_target = out_group`
            `vertices_source = out_group`, `vertices_target = in_group`

See:

Bell, JR. Subgroup centrality measures. Network Science. 2014;2(2):277-297. 
doi:10.1017/nws.2014.15



# Construct

```
subgroup_centrality_bells(
    graph::AbstractGraph;
    vertices_source::Union{Vector{Int64}, Nothing} = nothing,
    vertices_target::Union{Vector{Int64}, Nothing} = nothing,
)
```

##  Function Arguments

- `graph`: AbstractGraph 


##  Function Arguments

- `vertices_source`: set of vertex indices to treat as sources 
- `vertices_target`: set of vertex indices to treat as targets

"""
function subgroup_centrality_bells(
    graph::AbstractGraph;
    vertices_source::Union{Vector{Int64}, Nothing} = nothing,
    vertices_target::Union{Vector{Int64}, Nothing} = nothing,
)
    # get 
    vertices_source = isa(vertices_source, Nothing) ? collect(vertices(graph)) : vertices_source
    vertices_target = isa(vertices_target, Nothing) ? collect(vertices(graph)) : vertices_target
    vec_betweenness_centrality = zeros(Float64, nv(graph))

    for i in vertices_source
        b_i!(
            graph,
            i,
            vec_betweenness_centrality;
            vertices_target = vertices_target,
        )
    end

    !is_directed(graph) && (vec_betweenness_centrality ./= 2)
    
    return vec_betweenness_centrality
end





