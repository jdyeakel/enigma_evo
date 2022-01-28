using Distributed
using DataFrames
using Images
using UnicodePlots

@everywhere using AxisArrays
@everywhere using Combinatorics
@everywhere using LinearAlgebra
# @everywhere using Distributed
@everywhere using SharedArrays
@everywhere using SparseArrays
@everywhere using DataFrames
@everywhere using Distributions
@everywhere using Images
@everywhere using SpecialFunctions
@everywhere using LightGraphs
@everywhere using RCall
# @everywhere using HDF5
@everywhere using JLD2

if homedir() == "/home/z840"

    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/smartpath.jl")
    # @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/preamble_defs.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/intmatrixv4integer.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/intbool.jl")

    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/intmatrixv5.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/intadd.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/intfind.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/intfind_out.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/intfind_inout.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/potcol2.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/potextinct2.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/potsecextinct2.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/potobextinct2.jl")
    # @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/strength.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/strength2.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/secexteval.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/mutation.jl")


    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/assemblyeco.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/assemblyevo2.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/assemblyevo_diverse.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/assemblystate.jl")

    #Analysis Calculations
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/structure.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/dynstructure.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/sortassembly.jl")

    #Analysis functions
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/trophicalc2.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/roverlap.jl")

    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/nichemodelweb.jl")


else

    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/smartpath.jl")

    #Interaction matrix
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/intmatrixv4integer.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/intbool.jl")

    # @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/preamble_defs.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/intmatrixv5.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/intadd.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/intfind.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/intfind_out.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/intfind_inout.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/potcol2.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/potextinct2.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/potsecextinct2.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/potobextinct2.jl")
    # @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/strength.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/strength2.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/secexteval.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/mutation.jl")

    
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/assemblyeco.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/assemblyevo2.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/assemblyevo_diverse.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/assemblystate.jl")


    #Analysis Calculations
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/structure.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/dynstructure.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/sortassembly.jl")

    #Analysis functions
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/trophicalc2.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/roverlap.jl")

    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/nichemodelweb.jl")
end
