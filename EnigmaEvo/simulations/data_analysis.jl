simulation_name = "vary_cn_no_engineers_20_000_it_restrict_col";

avg_freqe, avg_freqe_pool, avg_freqn, avg_freqn_pool, avg_sprich, avg_pool  = 
    average_time_series(["freqe","freqe_pool","freqn","freqn_pool","sprich","pool"],
                        simulation_name,"cn",0:0.5:5,50)

jldsave("data/$(simulation_name)/averages.jld2",true;avg_freqe, avg_freqe_pool, avg_freqn, avg_freqn_pool, avg_sprich, avg_pool)

avg_freqe, avg_freqe_pool, avg_freqn, avg_freqn_pool, avg_sprich, avg_pool = 
    load("EnigmaEvo/data/$(simulation_name)/averages.jld2", "avg_freqe", "avg_freqe_pool",
        "avg_freqn", "avg_freqn_pool", "avg_sprich", "avg_pool");

itterations = 1:20_000

lineplot(itterations,hcat(avg_sprich[:,[1,end]],avg_pool[:,[1,end]]),xlabel="itteration",ylabel="species richness",title="Species richness of new model")
lineplot(itterations,hcat(avg_freqe[:,[1,end]],avg_freqe_pool[:,[1,end]]),xlabel="itteration",ylabel="freq eat",title="Frequency of eat interactions of new model")
lineplot(itterations,hcat(avg_freqn[:,[1,end]],avg_freqn_pool[:,[1,end]]),xlabel="itteration",ylabel="freq need",title="Frequency of need interactions of new model")

freqe, freqn,freqe_pool,freqn_pool, sprich,pool = load("EnigmaEvo/data/vary_cn_no_engineers/cn=5.0_repet=2.jld2", "freqe", "freqn","freqe_pool","freqn_pool", "sprich","pool");
itterations = 1:20000

lineplot(itterations,hcat(sprich,pool),xlabel="itteration",ylabel="species richness",title="Species richness of new model")
lineplot(itterations,hcat(freqe,freqe_pool),xlabel="itteration",ylabel="freq eat",title="Frequency of eat interactions of new model")
lineplot(itterations,hcat(freqn,freqn_pool),xlabel="itteration",ylabel="freq need",title="Frequency of need interactions of new model")


paramName = "nBasalRes"
simulationName = "vary_$(paramName)_50_000_it";        #specify the name of the simulation

specRich_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "species richness",legend = :topright);
meanEats_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "average amount of eat interactions",legend = :topright);
meanNeeds_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "average amount of need interactions",legend = :topright);

finalSpecRichPlt = plot(xlabel="number of basal resources",ylabel="mean spec richness in last 500 itterations", legend=nothing);

param_vals = 11:10:101
for nBasalRes in param_vals
    for rep in 1:10
        filename = "$(paramName)=$(nBasalRes)_repet=$(rep).jld2"
        sd = load("data/$(simulationName)/$(filename)", "simulationData");
        if  sd.clock[end] != 0. && maximum(sd.clock) < 1000
            plot!(specRich_plt,sd.clock,sd.specRich,color=RGB(nBasalRes/param_vals[end],0,0),label = "$(paramName) = $(nBasalRes)");
            plot!(meanEats_plt,sd.clock,sd.meanEats,color=RGB(nBasalRes/param_vals[end],0,0),label = "$(paramName) = $(nBasalRes)");
            plot!(meanNeeds_plt,sd.clock,sd.meanNeeds,color=RGB(nBasalRes/param_vals[end],0,0),label = "$(paramName) = $(nBasalRes)");
            scatter!(finalSpecRichPlt,[nBasalRes], [mean(sd.specRich[45_000:end])]);
        end
    end
end

Plots.savefig(specRich_plt,"data/$simulationName/plots/specRich_plt.png");
Plots.savefig(meanEats_plt,"data/$simulationName/plots/meanEats_plt.png");
Plots.savefig(meanNeeds_plt,"data/$simulationName/plots/meanNeeds_plt.png");
Plots.savefig(finalSpecRichPlt,"data/$simulationName/plots/finalSpecRichPlt.png");

display(specRich_plt)
display(meanNeeds_plt)
display(meanEats_plt)

display(finalSpecRichPlt)



node = getnode(phyloTree,getparent(phyloTree,"51v2"))
node.data["evolution"]


maxTrophLevels = maximum.(sd.trophLevels)
findmax(maxTrophLevels)

colnet_1950 = recreatecolnetdiverse(sd,1950)
colnet_1950.spec

eatMatrix = ENIgMaGraphs.convertToEatMatrixNonReduced(col_3770)
eatMatrixReduced = ENIgMaGraphs.convertToEatMatrix(col_3770)
R"""
    library(MASS)  
    library(NetIndices)
    rtl<-TrophInd(t($medRedEatMat))
"""
@rget rtl;
trophLevels = rtl[:,:TL]
maximum(trophLevels)

"$(sort(trophLevels))"

inds = sort(collect(union(col_3770.spec,col_3770.basalRes)))


medRedEatMat = eatMatrix[inds,inds]


paramName = "sqrt.(rPrimExt,rSecExt)"
simulationName = "varyExtinctionsForTrophLevel";        #specify the name of the simulation
sqrtParamVals = 1:.5:10;       #specify the parameter values that shall be simulated
repets = [1]        #specify the number of repetitions per parameter value
           #should the data be compressed before storing?
loop_vars = [(sqrtRPrimExt,sqrtRSecExt,repetition) for sqrtRPrimExt in sqrtParamVals for sqrtRSecExt in sqrtParamVals for repetition in repets];

using SharedArrays
numParams = length(sqrtParamVals)
heatMapMax = SharedArray{Float64}((numParams,numParams));
heatMapMean = SharedArray{Float64}((numParams,numParams));
@distributed for (sqrtRPrimExt,sqrtRSecExt,repetition) in loop_vars
    sd = load("EnigmaEvo/data/$(simulationName)/$(paramName)=($(sqrtRPrimExt),$(sqrtRSecExt))_repet=$repetition.jld2",  "simulationData")
    indSec = findfirst(==(sqrtRSecExt),sqrtParamVals)
    indPrim = findfirst(==(sqrtRPrimExt),sqrtParamVals)

end
using PlotlyJS
gr()
PlotlyJS.plot(PlotlyJS.heatmap(x = paramVals .^2, y = paramVals .^2, z=heatMapMax))
UnicodePlots.heatmap((heatMapMax))
UnicodePlots.heatmap((heatMapMean))
data = [1 20 30; 20 1 60; 30 60 1]
Plots.heatmap(heatMapMax)


data = [1 20 30; 20 1 60; 30 60 1]

PlotlyJS.plot(PlotlyJS.heatmap(z=data))



function analyse()
    paramName = "(rPrimExt,rSecExt)"
    simulationName = "varyExtinctionsForTrophLevel"; 
    compress = true
    paramVals = (1:.5:10).^2;       #specify the parameter values that shall be simulated
    numParams = length(paramVals)
    nRepets = 50
    repets = 1:nRepets
            
    loop_vars = [(primInd,rPrimExt,secInd,rSecExt,repetition) for (primInd,rPrimExt) in enumerate(paramVals)
        for (secInd,rSecExt) in enumerate(paramVals) for repetition in repets];

    monitoredParams = 
    specRichHeatMapMax = SharedArray{Float64}((numParams,numParams,nRepets));
    specRichHeatMapMean = SharedArray{Float64}((numParams,numParams,nRepets));
    @sync @distributed for (primInd,rPrimExt,secInd,rSecExt,repetition) in loop_vars

        sd = load("data/$(simulationName)/$(paramName)=($(rPrimExt),$(rSecExt))_repet=$repetition.jld2","simulationData")
        
        specRich = sd.specRich[9500:end]
        specRichHeatMapMax[primInd,secInd,repetition] = mean(maximum.(specRich))
        specRichHeatMapMean[primInd,secInd,repetition] = mean(mean.(specRich))
    end

    jldsave("data/$(simulationName)/results.jld2",compress; heatMapMax, heatMapMean)

    plotlyjs()
    maxPlot = Plots.heatmap(paramVals, paramVals, dropdims(mean(specRichHeatMapMax,dims=3),dims=3),
        size = (1280,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
        title = "Average maximal trophic level in itterations 9500 to 10000")
    meanPlot = Plots.heatmap(paramVals, paramVals, dropdims(mean(specRichHeatMapMean,dims=3),dims=3),
        size = (1280,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
        title = "Average mean trophic level in itterations 9500 to 10000")

    Plots.savefig(maxPlot,"data/$simulationName/plots/maxTrophLevelPlot.html");
    Plots.savefig(meanPlot,"data/$simulationName/plots/meanTrophLevelPlot.html");
end



paramName = "(rPrimExt,rExt)"
simulationName = "heatmapPrimExtVsGlobExt";        #specify the name of the simulation
compress::Bool = true;  #should the data be compressed before storing?


primVals = [.1,.5,1,1.5,2,3,5,10,15,20,40]       #specify the parameter values that shall be simulated
globVals = [.001,.005,.01,.016,.3,.35,.45,.7,1.,2.,5.]
numPrimVals = length(primVals)
numGlobVals = length(globVals)

nRepets = 50
repets = 1:nRepets


@enum NetworkGrowthType begin
    uncertain = 0
    boundedUncertain = 1
    boundedStable = 2
    boundedUnstable = 3
    unboundedUncertain =4
    undboundedLinear = 5
    unboundedExponential = 6
    trivial = 7
end


maxits = 10_000
#prepare everything for a simulation consisting of the variation of a parmeter
paramName = "(rPrimExt,rSecExt)"
simulationName = "varyExtinctionsForTrophLevel";        #specify the name of the simulation
mkpath("data/$(simulationName)/runs");    #make a folder with that name in the Data folder including plots subfolder
mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder including plots subfolder
compress::Bool = true;  #should the data be compressed before storing?


paramVals = (1:0.5:10).^2#0.5:.5:6.5       #specify the parameter values that shall be simulated
numParams = length(paramVals)
nRepets = 50
repets = 1:nRepets

repetition = 2

plots = Array{Plots.Plot}(undef,numParams,numParams)
nets = Array{ENIgMaSimulationData}(undef,numParams,numParams)
using SharedArrays
bounded = SharedArray{Bool}((numParams,numParams,nRepets))
netGrowthTypes = SharedArray{NetworkGrowthType}((numParams,numParams,nRepets))


const maxSlope = 0.1
const maxRecursionDepth = 5

using GLM
using DataFrames

@everywhere function checkBoundedness(sd,rates,bounded,primInd,rPrimExt,secInd,rSecExt,repet,recursionDepth=0;offset = 0)
    linReg = lm(@formula(specRich ~ clock), DataFrame(clock = sd.clock[offset:end], specRich = sd.specRich[offset:end]))

    coefTable = GLM.coeftable

    lower95 = coefTable[5][2]
    upper95 = coefTable[6][2]

    if lower95 > maxSlope
        bounded[primInd,secInd,parameter] = false
        netGrowthTypes[primInd,secInd,parameter] = unboundedUncertain
    elseif upper95 < maxSlope# && lower95 > -maxSlope
        bounded[primInd,secInd,parameter] = true
        netGrowthTypes[primInd,secInd,parameter] = boundedUncertain
    else
        if recursionDepth < maxRecursionDepth
            sd,_ = assemblyevo(poolnet, rates, maxits, cn,cn,ce,cpred, diverse, restrict_colonization, sd.colnet)
            checkBoundedness(sd,rates,bounded,primInd,rPrimExt,secInd,rSecExt,repet,recursionDepth+1)
        else
            netGrowthTypes[primInd,secInd,parameter] = uncertain
        end
    end
end

@sync @distributed for (primInd,rPrimExt) in enumerate(paramVals), (secInd,rSecExt) in enumerate(paramVals), repetition in repets
    sd,rates0 = load("data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rSecExt))_repet=$repetition.jld2", "simulationData","rates0")

    checkBoundedness(sd,rates,bounded,primInd,rPrimExt,secInd,rSecExt,repet;offset = 2_000)
    #nets[primInd,secInd] = sd
    #plots[primInd,secInd] = plot(sd.clock,sd.specRich,title="(rPrimExt,rSecExt) = $((rPrimExt,rSecExt))",ylabel = "specRich", xlabel = "clock", show = false)
end

boundedPlot = Plots.heatmap(string.(paramVals), string.(paramVals), dropdims(count(bounded,dims=3),dims=3),
size = (720,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
title = "Number of runs categorized as bounded.", xticks = :all,
yticks = :all, xrotation = 60)

Plots.savefig(boundedPlot,"data/$simulationName/plots/boundedPlot.html");



plots

summaryPlt = plot(plots..., size = (3500, 1500), legend = false)

Plots.savefig(summaryPlt,"data/$(simulationName)/plots/specRichTimeSeriesOvervies_2.html")

rPrimExt,rSecExt = 1.,1.

GLM.coeftable(linReg).cols[5][2]