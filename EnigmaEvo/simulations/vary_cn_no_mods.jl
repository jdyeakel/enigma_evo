if homedir() == "/home/z840"    #for downward compatibility ;)
    localpath::String = "$(homedir())/2019_Lego_Evo/EnigmaEvo/src/";
elseif isfile("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl")
    localpath::String = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/";
else                            #othserwise use relative paths
    localpath::String = "../src/";
end;

#using Revise    #helps with debugging in REPL (automatically tracks changes eg in files included with "includet" (include and track))
include(localpath*"loadfuncs.jl"); 

function simulation()   #supposedly its better to wrap stuff in functions and not use global variables if possible(not sure if typed ones are a problem)
    include(localpath*"set_up_params.jl");

    #prepare everything for a simulation consisting of the variation of a parmeter
    param_name = "cn"
    simulation_name = "vary_$(param_name)_no_engineers";        #specify the name of the simulation
    mkpath("data/$(simulation_name)");    #make a folder with that name in the Data folder

    param_vals = 1:10:51;        #specify the parameter values that shall be simulated
    repetitions_per_param = 5;      #specify the number of repetitions per parameter value

    compress::Bool = true;               #should the data be compressed before storing?
    loop_vars = collect((param,repetition) for param in param_vals for repetition in 1:repetitions_per_param);

    @distributed for (param,repetition) in loop_vars
        initpoolnet::ENIgMaGraph = setuppool(S,lambda,SSprobs,SOprobs,diverse);

        # EVOLUTIONARY VERSION
        sd = assemblyevo(initpoolnet, rates0, maxits, cn,cm,ce,cpred, diverse, restrict_colonization, logging);
        
        jldsave("data/$(simulation_name)/$(param_name)=$(param)_repet=$repetition.jld2",compress; simulationData = sd, rates0, maxits, param,cm,ce,cpred, diverse, restrict_colonization, logging,S,lambda,SSprobs,SOprobs)
    end
end

simulation()