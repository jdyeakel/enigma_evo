if homedir() == "/home/z840"    #for downward compatibility ;)
    localPath::String = "$(homedir())/2019_Lego_Evo/EnigmaEvo/src/";
elseif isfile("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl")
    localPath::String = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/";
else                            #othserwise use relative paths
    localPath::String = "../src/";
end;

#using Revise    #helps with debugging in REPL (automatically tracks changes eg in files included with "includet" (include and track))
include(localPath*"loadfuncs.jl");

function simulation()   #supposedly its better to wrap stuff in functions and not use global variables if possible(not sure if typed ones are a problem)
    include(localPath*"set_up_params.jl");
    maxits = 50_000
    #prepare everything for a simulation consisting of the variation of a parmeter
    paramName = "nBasalRes"
    simulationName = "vary_$(paramName)_50_000_it";        #specify the name of the simulation
    mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder
    compress::Bool = true;
    jldsave("data/$(simulationName)/parameters.jld2",compress;
        rates0, maxits, cm,ce,cpred, diverse, restrict_colonization,
        logging,S,lambda,SSprobs,SOprobs,nBasalRes)

    paramVals = 11:10:101;       #specify the parameter values that shall be simulated
    repetsPerParam = 10;        #specify the number of repetitions per parameter value
               #should the data be compressed before storing?
    loop_vars = collect((param,repetition) for param in paramVals for repetition in 1:repetsPerParam);

    @distributed for (nBasalRes,repetition) in loop_vars
        initpoolnet::ENIgMaGraph = setUpPool(S,lambda,nBasalRes,SSprobs,SOprobs,diverse);

        # EVOLUTIONARY VERSION
        (sd,_) = assemblyevo(initpoolnet, rates0, maxits, cm,cn,ce,cpred, diverse, restrict_colonization);
        
        jldsave("data/$(simulationName)/$(paramName)=$(nBasalRes)_repet=$repetition.jld2",compress; simulationData = sd)
    end
end

simulation()