using Distributed
#using DataFrames
#using Images
using UnicodePlots
using Plots

#@everywhere using AxisArrays
#@everywhere using Combinatorics
#@everywhere using LinearAlgebra
#@everywhere using Distributed
#@everywhere using SharedArrays
#@everywhere using SparseArrays
#@everywhere using DataFrames
@everywhere using OffsetArrays
@everywhere using Distributions
@everywhere using Images
#@everywhere using SpecialFunctions
#@everywhere using LightGraphs
@everywhere using RCall
@everywhere using Random
# @everywhere using HDF5
@everywhere using Graphs
@everywhere using JLD2
@everywhere using CodecZlib #used for compression by JLD2
@everywhere using Phylo

@everywhere include("ENIgMaGraphs.jl")
@everywhere using .ENIgMaGraphs
@everywhere include("setuppool.jl")
@everywhere include("assemblyevo3.jl")
@everywhere include("average_time_series.jl")
@everywhere include("extinction_size.jl")
include("plots.jl")

#=
@everywhere include("smartpath.jl")


#Interaction matrix
@everywhere include("intmatrixv4integer.jl")
@everywhere include("intbool.jl")

@everywhere include("preamble_defs.jl")
@everywhere include("intmatrixv5.jl")
@everywhere include("intadd.jl")
@everywhere include("intfind.jl")
@everywhere include("intfind_out.jl")
@everywhere include("intfind_inout.jl")
@everywhere include("potcol2.jl")
@everywhere include("potextinct2.jl")
@everywhere include("potsecextinct2.jl")
@everywhere include("potobextinct2.jl")
@everywhere include("strength.jl")
@everywhere include("strength2.jl")
@everywhere include("secexteval.jl")
@everywhere include("mutation.jl")


@everywhere include("assemblyeco.jl")
@everywhere include("assemblyevo2.jl")
@everywhere include("assemblyevo_diverse.jl") 
@everywhere include("assemblystate.jl")


#Analysis Calculations
@everywhere include("structure.jl")
@everywhere include("dynstructure.jl")
@everywhere include("sortassembly.jl")

#Analysis functions
@everywhere include("trophicalc2.jl")
@everywhere include("roverlap.jl")
=#