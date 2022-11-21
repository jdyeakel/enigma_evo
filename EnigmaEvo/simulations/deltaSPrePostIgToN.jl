if homedir() == "/home/z840"    #for downward compatibility ;)
    srcPath::String = "$(homedir())/2019_Lego_Evo/EnigmaEvo/src/";
elseif isfile("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl")
    srcPath::String = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/";
else                            #othserwise use relative paths
    srcPath::String = "../src/";
end;

using SharedArrays
using DataFrames
#using Revise    #helps with debugging in REPL (automatically tracks changes eg in files included with "includet" (include and track))
include(srcPath*"loadfuncs.jl"); 

function simulation(srcPath)   #supposedly its better to wrap stuff in functions and not use global variables if possible (with typed globals maybe only to save memory after distributed as globals are copied in each process and not deleted after completion)
    include(srcPath*"set_up_params.jl");

    #prepare everything for a simulation consisting of the variation of a parmeter
    simulation_name = "deltaSPrePostIgToN";        #specify the name of the simulation
    mkpath("data/$(simulation_name)");    #make a folder with that name in the Data folder

    maxits = 10_000
    repetitions = 1000;      #specify the number of repetitions per parameter value

    timeWindows = [.1,.25,.5,.8,1.,2.,5.,10.]
    offset = 2000

    compress::Bool = true;  j
    deltaSPrePostDicts = @distributed (vcat) for repetition in 1:repetitions
        
        initpoolnet = setuppool(S,lambda,SSprobs,SOprobs,diverse);

        # EVOLUTIONARY VERSION
        sd = assemblyevo(initpoolnet, rates0, maxits, cn,cm,ce,cpred, diverse, restrict_colonization, logging);
        
        jldsave("data/$(simulation_name)/repet=$repetition.jld2",compress; simulationData = sd, rates0, maxits,cm,ce,cpred, diverse, restrict_colonization, logging,S,lambda,SSprobs,SOprobs)
        
        evos = findall(ev -> isMutationType(ev,ignoreInteraction,needInteraction), sd.events)
        deltaSPrePostDict = Dict{Float64,Tuple{Array{Int},Array{Int}}}()
        for deltaT in timeWindows
             deltaSPrePostDict[deltaT] = deltaSPrePostEvents(sd.sprich, evos, sd.clock, deltaT, offset=offset )
        end

        deltaSPrePostDict
    end

    deltaSPrePostDict = Dict{Float64,Tuple{Array{Int},Array{Int}}}()
    for deltaT in timeWindows
        deltaTPre = vcat([ deltaTPre for (deltaTPre,_) in (deltaSPrePostDicts[repet][deltaT] for repet in 1:repetitions) ]...)
        deltaTPost = vcat([ deltaTPre for (_,deltaTPost) in (deltaSPrePostDicts[repet][deltaT] for repet in 1:repetitions) ]...)
        deltaSPrePostDict[deltaT] = (deltaTPre,deltaTPost)
    end

    jldsave("data/$(simulation_name)/results.jld2",compress;deltaSPrePostDict)
end

simulation(srcPath)

deltaSPrePostDict = load("data/deltaSPrePostIgToN/results.jld2","deltaSPrePostDict")