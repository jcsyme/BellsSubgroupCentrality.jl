###   ORDERED LOADING

module BellsSubgroupCentrality

using Graphs

# export key functions\
export get_subgroup_centrality,
       subgroup_centrality_bells


# read from file
dir_load = @__DIR__
include(joinpath(dir_load, "bells_centrality.jl"))

end