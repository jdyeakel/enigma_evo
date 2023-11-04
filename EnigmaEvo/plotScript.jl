#script for the final plots of the report

if homedir() == "/home/z840"    #for downward compatibility ;)
    localPath::String = "$(homedir())/2019_Lego_Evo/EnigmaEvo/src/";
elseif isfile("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl")
    localPath::String = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/";
else                            #othserwise use relative paths
    localPath::String = "src/";
end;

Base.InvasiveLinkedList = Base.IntrusiveLinkedList  #this line works around a renaming issue caused by an update of julia libraries

#load librarys
include(localPath*"loadfuncs.jl");

using SharedArrays
@everywhere import NaNStatistics

NO_STANDALONE = true;   #this variable is used to let the plot files know they are called by a script and thus don't need to load librarys
#load plot functions
include("simulations/setUpPlotProperties.jl")
include("simulations/specRichAndTrophLvlOverPrimExt.jl")
include("simulations/rEvoLinePlot.jl")
include("simulations/similarityAnalysis.jl")

#create plots
legendPosition = :rt
rPrimExtLinePlot(true,true,plotLib = :Makie)
legendPosition = :lt
rEvoLinePlot(true,true,plotLib = :Makie)
doSimilarityAnalysis(true,true,plotLib = :Makie)