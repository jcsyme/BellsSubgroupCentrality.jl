
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
Slow implementation of betweenness centrality that can be used to implment Bell's centrality.

Specify groups using vertices_source and vertices_target
"""
function betweenness_centrality_bells(
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



