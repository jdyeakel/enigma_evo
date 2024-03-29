#some how tos at end of file
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

#Random.seed!(2456526);
#create random pool network
poolnet::ENIgMaGraph = setUpPool(S,lambda,nBasalRes,SSprobs,SOprobs,diverse);

# run a simulation with parameters given (always use a freshly initialized poolnet as the poolnet is changed during assembly)
#results are stored in a ENIgMaSimulationData subtype (sd shorthand alias) (see assemblyevo3.jl for definitions)
# a second return value can be used to return extra data, that might be temporarily be of interest 
@time simulationData,_ = sd,extraData = 
    assemblyevo(poolnet, rates0, maxits, cm,cn,ce,cpred, diverse, restrict_colonization, createLog = true);

#plot some results:
plotlyjs();    #use the plotlyjs backend for interactivity
#gr() use the gr backend for faster plots

#Plot some time series and the distribution of extinction sizes ignoring the first 2000 itterations
plotSimulation(sd,offset=2000,show=true)

#plot the phylogeny
plotPhylogeny(sd.phyloTree,sorted=true)

#plot mean and max trophic level (as computation via R is slow the trophic level is not computed every itteration and sometimes only in the last few thousand itterations at all)
plot(sd.clock[28_500:30:end],[mean.(sd.trophLevels),maximum.(sd.trophLevels)], xlabel = "clock time", ylabel="trophic level", label=["mean trophic level" "maximal trophic level"])

#plot some more properties
plot(sd.clock,[sd.specRich,sd.pool,sd.nColonizers,sd.specRich + sd.nColonizers,sd.nSecExtSpec,sd.nPrimExtSpec, sd.specRich - sd.nSecExtSpec - sd.nPrimExtSpec], size = (1920,1080),
    label = ["species richness" "pool spec richness" "#potential colonizers" "specRich + colonizers"  "#secondary ext. species" "#primary ext. species" "specRich - nPrimExt - n SecExt"])

#following paragraph: early attempts to track triggers of extinction cascades; can probably be ignored as not very well documented 
#what happens after a certain sort of events to species richness?
#here: what happens to the species richenss S after ignore to need mutations (here within .8 time units)
evos = findall(ev -> isMutationType(ev,ignoreInteraction,needInteraction), sd.events)
ΔSPre,ΔSPost = deltaSPrePostEvents(sd.specRich, evos, sd.clock, .8, offset=1000 )
Plots.histogram(ΔSPre, alpha=.5, normalize=true,yaxis=:log,xaxis="ΔS pre and post event", label="pre event", title = "Ignore to need mutations")
Plots.histogram!(ΔSPost, alpha=.5, normalize=true,yaxis=:log, label="post event")

# save everything you want to use later in file 
compress = true;    #should data be compressed?
#put all variable to be saved after semicolon, order doesnt matter, name of variable will be identifier in file
jldsave("tutorial.jld2",compress;simulationData,rates0)

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

#reconstruct the state of the colony at a certain itteration (not finally tested)
colnet_it_100 = recreatecolnetdiverse(sd,100)

#get number of species in a network
numspec(colnet_it_100)

#for some fuctions there is : more information available press ? in REPL(console), then enter function name (or hover over in vs code)
# try: ? recreatecolnetdiverse

#outdated # average one or multiple time series from a simulation
#avg_sprich =  average_time_series("sprich","vary_cn_no_engineers","cn",0:.1:5,2)
#avg_sprich, avg_freqn =  average_time_series(["freqn","freqe"],"vary_cn_no_engineers","cn",0:.1:5,2)
