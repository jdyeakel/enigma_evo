using Distributed
using DataFrames
using Images

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

    #Interaction matrix
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/intmatrixv3.jl")

		#Alternative 4x4 Interaction matrix
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/intmatrixv4.jl")

    #Community dynamics
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/preamble_defs.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/assemblyeco.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/assemblyevo.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/assemblystate.jl")

    #Analysis Calculations
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/structure.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/dynstructure.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/sortassembly.jl")

    #Analysis functions
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/trophicalc2.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/roverlap.jl")
    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/potcol.jl")

    @everywhere include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/nichemodelweb.jl")


else

    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/smartpath.jl")

    #Interaction matrix
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/intmatrixv3.jl")

		#Alternative 4x4 Interaction matrix
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/intmatrixv4.jl")

    #Community dynamics
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/preamble_defs.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/assemblyeco.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/assemblyevo.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/assemblystate.jl")


    #Analysis Calculations
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/structure.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/dynstructure.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/sortassembly.jl")

    #Analysis functions
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/trophicalc2.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/roverlap.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/potcol.jl")

    @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/nichemodelweb.jl")
end
