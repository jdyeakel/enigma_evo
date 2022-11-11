# get path to local src directory
if homedir() == "/home/z840"    #for downward compatibility ;)
    localpath::String = "$(homedir())/2019_Lego_Evo/EnigmaEvo/src/";
elseif isfile("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl")
    localpath::String = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/";
else                            #othserwise use relative paths
    localpath::String = "src/";
end;

#load all necessary functions
include(localpath*"loadfuncs.jl"); 

#setup all parameters
include(localpath*"set_up_params.jl");

#create random pool network
poolnet::ENIgMaGraph = setuppool(S,lambda,SSprobs,SOprobs);

# run a simulation with parameters given (always use a freshly initialized poolnet as the poolnet is changed during assembly)
@time poolnet,colnet,phyloTree,sprich,rich,pool,mstrength,evolvedstrength,clock,CID,maxids,glob_ext_spec,mutstep,freqe,freqn,freqe_pool,freqn_pool,events =
    simulation_data = assemblyevo(poolnet, rates0, maxits, cn,cn,ce,cpred, diverse, restrict_colonization, logging);

# save everything you want to use later in file 
plotcompress = true;    #should data be compressed?
#put all variable to be saved after semicolon, order doesnt matter
jldsave("tutorial.jld2",compress;poolnet,colnet,sprich,clock,CID,glob_ext_spec,maxids)

#for loading use eg the following
poolnet_loaded, colnet_loaded = load("tutorial.jld2", "poolnet","colnet");

#How to:-------------------------------------------------------------------------------
#get species with id from network
spec = poolnet[2]
#get a species interactions (eat, feed, make or need)
needs = spec.need;
#check if a vertex(modifier or species) with id is in the network
exists = colnet.hasv[3]
#check if a species with id is in the network/if the vertex with that id is a species
isspec = colnet.hasspec[3]
#get number of species in a network
numspec(colnet_it_100)
#reconstruct the state of the colony at a certain itteration
colnet_it_100 = recreatecolnetdiverse(poolnet::ENIgMaGraph,100,CID[:,100],maxids[100],glob_ext_spec)
#for some fuctions there is : more information available press ? in REPL(console), then enter function name (or hover over in vs code)
# try: ? recreatecolnetdiverse

# average one or multiple time series from a simulation
avg_sprich =  average_time_series("sprich","vary_cn_no_engineers","cn",0:.1:5,2)
avg_sprich, avg_freqn =  average_time_series(["freqn","freqe"],"vary_cn_no_engineers","cn",0:.1:5,2)
