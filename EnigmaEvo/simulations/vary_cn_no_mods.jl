if homedir() == "/home/z840"    #for downward compatibility ;)
    localpath::String = "$(homedir())/2019_Lego_Evo/EnigmaEvo/src/";
elseif isfile("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl")
    localpath::String = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/";
else                            #othserwise use relative paths
    localpath::String = "../src/";
end;

#using Revise    #helps with debugging in REPL (automatically tracks changes eg in files included with "includet" (include and track))
include(localpath*"loadfuncs.jl"); 

include(localpath*"set_up_params.jl");

#prepare everything for a simulation consisting of the variation of a parmeter
simulation_name = "vary_cn_no_engineers";        #specify the name of the simulation
mkpath("data/"*simulation_name);    #make a folder with that name in the Data folder

param_vals = 3.8:0.1:5;        #specify the parameter values that shall be simulated
repetitions_per_param = 2;      #specify the number of repetitions per parameter value

compress::Bool = true;               #should the data be compressed before storing?
loop_vars = collect((cn,repetition) for cn in param_vals for repetition in 1:repetitions_per_param);

Threads.@threads for (cn,repetition) in loop_vars
    println(Threads.threadid())
    initpoolnet::ENIgMaGraph = setuppool(S,lambda,SSprobs,SOprobs);

    # EVOLUTIONARY VERSION
    poolnet,colnet,sprich,rich,pool,mstrength,evolvedstrength,clock,CID,maxids,glob_ext_spec,mutstep,freqe,freqn,events =
        assemblyevo(initpoolnet, rates0, maxits, cn,cm,ce,cpred, diverse, restrict_colonization, logging);
    
    jldsave("data/"*simulation_name*"/cn=$(cn)_repet=$repetition.jld2",compress;poolnet,colnet,sprich,rich,pool,mstrength,evolvedstrength,clock,CID,mutstep,freqe,freqn,events, rates0, maxits, cn,cm,ce,cpred, diverse, restrict_colonization, logging,S,lambda,SSprobs,SOprobs)
end