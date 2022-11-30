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
    maxits = 10_000
    #prepare everything for a simulation consisting of the variation of a parmeter
    paramName = "sqrt.(rPrimExt,rSecExt)"
    simulationName = "varyExtinctionsForTrophLevel";        #specify the name of the simulation
    mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder
    compress::Bool = true;
    jldsave("data/$(simulationName)/parameters.jld2",compress;
        rates0, maxits, cm,ce,cpred, diverse, restrict_colonization,
        logging,S,lambda,SSprobs,SOprobs,nBasalRes)

    sqrtParamVals = 1:.5:10;       #specify the parameter values that shall be simulated
    repets = [1]        #specify the number of repetitions per parameter value
               #should the data be compressed before storing?
    loop_vars = [(sqrtRPrimExt,sqrtRSecExt,repetition) for sqrtRPrimExt in sqrtParamVals for sqrtRSecExt in sqrtParamVals for repetition in repets];

    @distributed for (sqrtRPrimExt,sqrtRSecExt,repetition) in loop_vars
        initpoolnet::ENIgMaGraph = setUpPool(S,lambda,nBasalRes,SSprobs,SOprobs,diverse);
        rates0 = (rc = 1., rprimext = sqrtRPrimExt^2, rsecext = sqrtRSecExt^2, reo = 1., revo = 0.035, rext = .016);
        # EVOLUTIONARY VERSION
        (sd,_) = assemblyevo(initpoolnet, rates0, maxits, cn,cm,ce,cpred, diverse, restrict_colonization);
        
        jldsave("data/$(simulationName)/$(paramName)=$(rPrimExt,rSecExt)_repet=$repetition.jld2",compress; simulationData = sd)
    end
end

simulation()