if homedir() == "/home/z840"    #for downward compatibility ;)
    localPath::String = "$(homedir())/2019_Lego_Evo/EnigmaEvo/src/";
elseif isfile("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl")
    localPath::String = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/";
else                            #othserwise use relative paths
    localPath::String = "src/";
end;

#using Distributed
#addprocs(80)

include(localPath*"loadfuncs.jl");

using SharedArrays
@everywhere import NaNStatistics
import StatsPlots

function simulation(onlyPlots=false,fromResults=false;couplingFactor=.5)   #supposedly its better to wrap stuff in functions and not use global variables if possible(not sure if typed ones are a problem)
    if fromResults && !onlyPlots
        error("If only results file should be used it implies that only plots can be created.")
    end

    include("src/set_up_params.jl");
    maxits = 30_000
    #prepare everything for a simulation consisting of the variation of a parmeter
    paramName = "rPrimExt,rSecExt,evoInd"
    simulationName = "rEvoLinePlotHighGlobExt";        #specify the name of the simulation
    mkpath("data/$(simulationName)/runs");    #make a folder with that name in the Data folder including plots subfolder
    mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder including plots subfolder
    compress::Bool = true;  #should the data be compressed before storing?


    primVals = [3.5]       #specify the parameter values that shall be simulated
    nPrimVals = length(primVals)
    secVals = [1.,10.]
    nSecVals = length(secVals)
    evoVals = exp10.(range(log10(0.005), stop=log10(2), length=40))[1:30]
    nEvoVals = length(evoVals)
    evoInds = [6,28]#,37];    #hand picked indices of values of special interest 
    nEvoInds = length(evoInds)
    nRepets = 50
    repets = 1:nRepets

    jldsave("data/$(simulationName)/parameters.jld2",compress;
    rates0, maxits, cm,ce,cpred, diverse, restrict_colonization,
    logging,S,lambda,SSprobs,SOprobs,nBasalRes,primVals,nRepets, secVals,evoVals,couplingFactor)

    endVectorialZParams = [:trophLevels]
    allTimeScalarZParams = Symbol[:specRich, :pool, :meanEats, :meanNeeds,  :nPrimExtSpec, :nSecExtSpec, :nColonizers]
    zParams = [endVectorialZParams; allTimeScalarZParams]
    zParamLongNames = Dict{Symbol,String}(
        :trophLevels => "trophic level",
        :specRich => "species richness",
        :pool => "species richness in the pool",
        :meanEats => "mean number of eats per species",
        :meanNeeds => "mean number of needs per species",
        :nPrimExtSpec => "number of primary extinction candidates",
        :nSecExtSpec => "number of secondary extinction candidates",
        :nColonizers => "number of potential Colonizers"
    )

    if !fromResults
        maxResults = Dict{Symbol, SharedArray{Float64}}()
        meanResults = Dict{Symbol,SharedArray{Float64}}()
        for vecZParam in endVectorialZParams
            maxResults[vecZParam] = SharedArray{Float64}((nPrimVals,nEvoVals,nSecVals,nRepets))
            maxResults[vecZParam][:] .= NaN
        end
        for zParam in zParams
            meanResults[zParam] = SharedArray{Float64}((nPrimVals,nEvoVals,nSecVals,nRepets))
            meanResults[zParam][:] .= NaN
        end

        runFinished = SharedArray{Bool}((nPrimVals,nEvoVals,nSecVals,nRepets))

        loop_vars = [(primInd,rPrimExt,secInd,rSecExt,evoInd,rEvo,repetition) for (primInd,rPrimExt) in enumerate(primVals)
            for (secInd,rSecExt) in enumerate(secVals) for (evoInd,rEvo) in enumerate(evoVals) for repetition in repets];

        @sync @distributed for (primInd,rPrimExt,secInd,rSecExt,evoInd,rEvo,repetition) in loop_vars
            println("$((primInd,rPrimExt,secInd,rSecExt,evoInd,rEvo,repetition))")
            fileName = "data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rSecExt),$(evoInd))_repet=$repetition.jld2" 
            if onlyPlots
                if isfile(fileName)
                    sd = load(fileName, "simulationData")
                else
                    runFinished[primInd,evoInd,secInd,repetition] = false  #maybe redundant as SharedArray defaults to zeo or here false but not sure if that behaviour is guaranteed
                    continue
                end
            else
                rates0 = (rc = 1., rprimext = rPrimExt, rsecext = rSecExt, reo = 1., revo = rEvo, rext = rEvo*couplingFactor);
                initpoolnet::ENIgMaGraph = setUpPool(S,lambda,nBasalRes,SSprobs,SOprobs,diverse);
                (sd,_) = assemblyevo(initpoolnet, rates0, maxits, cm,cn,ce,cpred, diverse, restrict_colonization);

                jldsave(fileName,compress; simulationData = sd, rates0)
            end
            if sd.runFinished    #did the run end prematurely?
                runFinished[primInd,evoInd,secInd,repetition] = true
                for endVectorialZParamSymb in endVectorialZParams
                    endVectorialZParam = getfield(sd,endVectorialZParamSymb)
                    if any(isempty,endVectorialZParam)
                        maxResults[endVectorialZParamSymb][primInd,evoInd,secInd,repetition] = NaN
                        meanResults[endVectorialZParamSymb][primInd,evoInd,secInd,repetition] = NaN
                    else
                        maxResults[endVectorialZParamSymb][primInd,evoInd,secInd,repetition] = mean(maximum.(endVectorialZParam))
                        meanResults[endVectorialZParamSymb][primInd,evoInd,secInd,repetition] = mean(mean.(endVectorialZParam))
                    end
                end
                for allTimeScalarZParamSymb in allTimeScalarZParams
                    allTimeScalarZParam = getfield(sd,allTimeScalarZParamSymb)[(maxits-1500):end]
                    meanResults[allTimeScalarZParamSymb][primInd,evoInd,secInd,repetition] = mean(allTimeScalarZParam)
                end
            else
                runFinished[primInd,evoInd,secInd,repetition] = false
            end
        end
        jldsave("data/$(simulationName)/results.jld2",compress; maxResults, meanResults,runFinished)
    else
        maxResults, meanResults, runFinished = 
            load("data/$(simulationName)/results.jld2", "maxResults", "meanResults", "runFinished")
    end
    
    plotlyjs()
    xLabel = "evolutionary rate ($(1/couplingFactor)*global extinction rate)"
    seriesLabels = reshape(["rSecExt = $rSecExt" for rSecExt in secVals], (1,nSecVals))
    boxLabels = vec(repeat(reshape(["rEvo = $(evoVals[evoInd])" for evoInd in evoInds],(1,nEvoInds)),nRepets));
    for (primInd,rPrimExt) in enumerate(primVals)
        plotPath = "data/$simulationName/plots/rPrimExt=$(rPrimExt)"
        mkpath(plotPath)
        boxPlotPath = "data/$simulationName/plots/boxplots/logY/rPrimExt=$(rPrimExt)"
        mkpath(boxPlotPath)

        nFinishedRunsPlot = scatter(evoVals,  dropdims(sum(selectdim(runFinished,1,primInd),dims=3),dims=3),
            xlabel = xLabel, ylabel = "number of successfull runs used", legend = :bottomleft,
            xaxis = :log, title = "Number of successfull runs used", labels = seriesLabels)

        Plots.savefig(nFinishedRunsPlot,"$(plotPath)/nFinishedRunsPlot.svg");

        for vecZParam in endVectorialZParams
            maxPlot = plot(evoVals, dropdims(mean(selectdim(maxResults[vecZParam],1,primInd),dims=3),dims=3),
                xlabel = xLabel, ylabel = zParamLongNames[vecZParam],
                #title = "Average maximal $(zParamLongNames[vecZParam]) in itterations $(maxits - 1500) to $maxits",
                xaxis = :log, labels = seriesLabels, legend = :topleft)

            Plots.savefig(maxPlot,"$(plotPath)/max_$(String(vecZParam))Plot.svg");

            if any(isnan,maxResults[vecZParam])
                maxPlot = plot(evoVals, NaNStatistics.nanmean(selectdim(maxResults[vecZParam],1,primInd),dim=3),
                xlabel = xLabel, ylabel = zParamLongNames[vecZParam],
                #title = "Average maximal $(zParamLongNames[vecZParam]) in itterations $(maxits - 1500) to $maxits - ignoring NaNs.", 
                xaxis = :log, labels = seriesLabels, legend = :topleft)

                Plots.savefig(maxPlot,"$(plotPath)/max_$(String(vecZParam))Plot_ignoreNaNs.svg");
            end


            StatsPlots.violin(boxLabels, vec(transpose(maxResults[vecZParam][primInd,evoInds,1,:])),
                ylabel = zParamLongNames[vecZParam], side = :left, label = "rSecExt = $(secVals[1])", yaxis = :log )
            violinPlot = StatsPlots.violin!(boxLabels, vec(transpose(maxResults[vecZParam][primInd,evoInds,2,:])),
                ylabel = zParamLongNames[vecZParam], side = :right, label = "rSecExt = $(secVals[2])", yaxis = :log )
            Plots.savefig(violinPlot,"$boxPlotPath/max_$(vecZParam)ViolinPlot.svg")
        end

        for zParam in zParams
            meanPlot = plot(evoVals, dropdims(mean(selectdim(meanResults[zParam],1,primInd),dims=3),dims=3),
                xlabel = xLabel, ylabel = zParamLongNames[zParam],
                #title = "Average mean $(zParamLongNames[zParam]) in itterations $(maxits - 1500) to $maxits",
                xaxis = :log, yaxis = :log, labels = seriesLabels, legend = :topleft)

            Plots.savefig(meanPlot,"$(plotPath)/mean_$(String(zParam))Plot.svg");

            if any(isnan,meanResults[zParam])
                meanPlot = plot(evoVals, NaNStatistics.nanmean(selectdim(meanResults[zParam],1,primInd),dim=3),
                xlabel = xLabel, ylabel = zParamLongNames[zParam],
                #title = "Average mean $(zParamLongNames[zParam]) in itterations $(maxits - 1500) to $maxits - ignoring NaNs",
                xaxis = :log, yaxis = :log, labels = seriesLabels, legend = :topleft)

                Plots.savefig(meanPlot,"$(plotPath)/mean_$(String(zParam))Plot_ignoreNaNs.svg");
            end

            StatsPlots.violin(boxLabels, vec(transpose(meanResults[zParam][primInd,evoInds,1,:])),
                ylabel = zParamLongNames[zParam], side = :left, label = "rSecExt = $(secVals[1])", yaxis = :log )
            violinPlot = StatsPlots.violin!(boxLabels, vec(transpose(meanResults[zParam][primInd,evoInds,2,:])),
                ylabel = zParamLongNames[zParam], side = :right, label = "rSecExt = $(secVals[2])", yaxis = :log )
            Plots.savefig(violinPlot,"$boxPlotPath/mean_$(zParam)ViolinPlot.svg")

        end
        StatsPlots.violin(boxLabels, vec(transpose(meanResults[:specRich][primInd,evoInds,1,:] ./ meanResults[:pool][primInd,evoInds,1,:])),
            ylabel = "specRich / pool", side = :left, label = "rSecExt = $(secVals[1])" )
        violinPlot = StatsPlots.violin!(boxLabels, vec(transpose(meanResults[:specRich][primInd,evoInds,2,:] ./ meanResults[:pool][primInd,evoInds,2,:])),
            side = :right, label = "rSecExt = $(secVals[2])" )
        Plots.savefig(violinPlot,"$boxPlotPath/colonyPoolRatioViolinPlot.svg")
    end
end

#=
paramName = "rPrimExt,rSecExt,evoInd";
simulationName = "rEvoInDepth";


primVals = [3.5]
nPrimVals = length(primVals)
secVals = [1.,10.]
nSecVals = length(secVals)
evoVals = exp10.(range(log10(0.005), stop=log10(2), length=40));
evoInds = [6,28,37];
nEvoInds = length(evoInds)
nRepets = 50;


repets = 1:nRepets;



loop_vars = [(primInd,rPrimExt,secInd,rSecExt,evoIndInd,evoInd,repetition) for (primInd,rPrimExt) in enumerate(primVals)
    for (secInd,rSecExt) in enumerate(secVals) for (evoIndInd,evoInd) in enumerate(evoInds) for repetition in repets];



sds = Array{ENIgMaSimulationData}(undef,nPrimVals,nSecVals,nEvoInds,nRepets)
for (primInd,rPrimExt,secInd,rSecExt,evoIndInd,evoInd,repetition) in loop_vars
    fileName = "data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rSecExt),$(evoInd))_repet=$repetition.jld2" 
    if isfile(fileName)
        sd = load(fileName, "simulationData")
    else
        continue
    end
    sds[primInd,secInd,evoIndInd,repetition] = sd
end

jldsave("data/$(simulationName)/detailedResults.jld2",true; simulationDatasets = sds)

sds = load("data/$(simulationName)/detailedResults.jld2", "simulationDatasets");

linkType = :eat
primInd = 1
evoIndInd = 3
degreeDistPlots = Array{Plots.Plot}(undef,7,7)
for repet in 1:49
    StatsPlots.violin(["rEvo = $(evoVals[evoInds[evoIndInd]])"], ENIgMaGraphs.getDegrees(sds[primInd,1,evoIndInd,repet].colnet,linkType),
        ylabel = "Links of type $linkType degree", side = :left, label = "rSecExt = $(secVals[1])" )
    violinPlot = StatsPlots.violin!(["rEvo = $(evoVals[evoInds[evoIndInd]])"], ENIgMaGraphs.getDegrees(sds[primInd,2,evoIndInd,repet].colnet,linkType),
        ylabel = "Links of type $linkType degree", side = :right, label = "rSecExt = $(secVals[2])" )
    degreeDistPlots[repet] = violinPlot
end

for repet in 1:49
    #bar(ENIgMaGraphs.getDegreeDist(sds[primInd,1,evoIndInd,repet].colnet,linkType)...,xticks = :all,
    #    xlabel = "Links of type $linkType degree",ylabel = "count", alpha = .5, label = "rSecExt = $(secVals[1])" )
    degreePlot = bar(ENIgMaGraphs.getDegreeDist(sds[primInd,2,evoIndInd,repet].colnet,linkType)...,xticks = :all, #yaxis = :log,
        xlabel = "Links of type $linkType degree",ylabel = "count", alpha = .5, label = "rSecExt = $(secVals[2])")
    degreeDistPlots[repet] = degreePlot
end

summaryPlt = plot(degreeDistPlots..., size = (1920, 1080), legend = false)

endVectorialZParams = [:trophLevels]
allTimeScalarZParams = Symbol[:specRich, :pool, :meanEats, :meanNeeds,  :nPrimExtSpec, :nSecExtSpec, :nColonizers]
zParams = [endVectorialZParams; allTimeScalarZParams]
zParamLongNames = Dict{Symbol,String}(
    :trophLevels => "trophic level",
    :specRich => "species richness",
    :pool => "species richness in the pool",
    :meanEats => "mean number of eats per species",
    :meanNeeds => "mean number of needs per species",
    :nPrimExtSpec => "number of primary extinction candidates",
    :nSecExtSpec => "number of secondary extinction candidates",
    :nColonizers => "number of potential Colonizers"
)

maxResults, meanResults, runFinished = 
    load("data/$(simulationName)/results.jld2", "maxResults", "meanResults", "runFinished")

boxLabels = vec(repeat(transpose(["rEvo = $(evoVals[evoInd])" for evoInd in evoInds]),nRepets));
for (primInd,rPrimExt) in enumerate primVals
    boxPlotPath = "data/$simulationName/plots/boxplots/rPrimExt=$(rPrimExt)"
    mkpath(boxPlotPath)
    for zParam in zParams
        #maxResults[zParam] = maxResults[zParam][:,:,evoInds,]
        StatsPlots.violin(boxLabels, vec(transpose(meanResults[primInd,1,evoInds,:])),
            ylabel = zParamLongNames[zParam], side = :left, label = "rSecExt = $(secVals[1])" )
        violinPlot = StatsPlots.violin!(boxLabels, vec(transpose(meanResults[primInd,2,evoInds,:])),
            ylabel = zParamLongNames[zParam], side = :right, label = "rSecExt = $(secVals[2])" )
        savefig(violinPlot,"$boxPlotPath/mean_$(zParam)ViolinPlot.svg")
    end
end
=#