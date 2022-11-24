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
poolnet::ENIgMaGraph = setUpPool(S,lambda,nBasalResources,SSprobs,SOprobs,diverse);

# run a simulation with parameters given (always use a freshly initialized poolnet as the poolnet is changed during assembly)
@time simulationData = sd = #results are stored in a ENIgMaSimulationData subtype (sd shorthand alias)
    simulation_data = assemblyevo(poolnet, rates0, maxits, cn,cn,ce,cpred, diverse, restrict_colonization, createLog = logging);


#plot some results:
plotlyjs()    #use the plotlyjs backend for interactivity
#gr() use the gr backend for faster plots

#Plot some time series and the distribution of extinction sizes ignoring the first 2000 itterations
plot_simulation(sd,offset=2000,show=true)

#plot the phylogeny
plotPhylogeny(sd.phyloTree,sorted=true)


#use functions of the Graphs package on ENIgMaGraphs
#@enter strongly_connected_components(colnet)

#what happens after a certain sort of events to species richness?
#here: what happens after ignore to need mutations
evos = findall(ev -> isMutationType(ev,ignoreInteraction,needInteraction), sd.events)
ΔSPre,ΔSPost = deltaSPrePostEvents(ds.sprich, evos, ds.clock, .8, offset=1000 )
Plots.histogram(ΔSPre, alpha=.5, normalize=true,yaxis=:log,xaxis="ΔS pre and post event", label="pre event", title = "Ignore to need mutations")
Plots.histogram!(ΔSPost, alpha=.5, normalize=true,yaxis=:log, label="post event")
#bar(pairs(preDist), label = "pre event")
#bar!( pairs(postDist), label = "post event", alpha = .5)

# save everything you want to use later in file 
compress = true;    #should data be compressed?
#put all variable to be saved after semicolon, order doesnt matter, name of variable will be identifier in file
jldsave("tutorial.jld2",compress;simulationData)

#for loading use eg the following
simulationData_loaded = load("tutorial.jld2", "simulationData");

#continue old run by setting final colony as initial colony of new run
@time simulationData_continued = sd = #results are stored in a ENIgMaSimulationData subtype (sd shorthand alias)
    simulation_data = assemblyevo(simulationData_loaded.poolnet, rates0, maxits, cn,cn,ce,cpred, diverse, restrict_colonization, simulationData_loaded.colnet; createLog = logging);

#How to:-------------------------------------------------------------------------------
#get species with id from network
spec = sd.poolnet[2]
#get a species interactions (eat, feed, make or need)
needs = spec.need;
#check if a vertex(modifier or species) with id is in the network
exists = sd.colnet.hasv[3]
#check if a species with id is in the network/if the vertex with that id is a species
isspec = sd.colnet.hasspec[3]

#reconstruct the state of the colony at a certain itteration
colnet_it_100 = recreatecolnetdiverse(sd,100)

#get number of species in a network
numspec(colnet_it_100)

#for some fuctions there is : more information available press ? in REPL(console), then enter function name (or hover over in vs code)
# try: ? recreatecolnetdiverse

#outdated # average one or multiple time series from a simulation
#avg_sprich =  average_time_series("sprich","vary_cn_no_engineers","cn",0:.1:5,2)
#avg_sprich, avg_freqn =  average_time_series(["freqn","freqe"],"vary_cn_no_engineers","cn",0:.1:5,2)
