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

    categories = Symbol[:SDsPool,:SDsCol,:SDsColInPool,:SDsPoolWithoutCol,:SDsPoolOnlyEats,:SDsColOnlyEats,:SDsColInPoolOnlyEats,:SDsPoolWithoutColOnlyEats]
    categoryLongNames = Dict{Symbol,String}(
        :SDsPool = "Species in pool",
        :SDsCol = "Species in colony",
        :SDsColInPool = "Species in colony as subset of pool",
        :SDsPoolWithoutCol = "Species in pool without species in colony",
        :SDsPoolOnlyEats = "Species in pool reduced to eats",
        :SDsColOnlyEats = "Species in colony reduced to eats",
        :SDsColInPoolOnlyEats = "Species in colony as subset of pool reduced to eats",
        :SDsPoolWithoutColOnlyEats = "Species in pool without species in colony reduced to eats"
    )

    if !fromResults
        SDmeans = Dict{Symbol, SharedArray{Float64}}()                                                                                                                                            bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
        SDstds = Dict{Symbol, SharedArray{Float64}}()                                                                                                                                            bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
        for category in categories
            SDmeans[category] = SharedArray{Float64}((nPrimVals,nEvoVals,nSecVals,nRepets))
            SDstds[category] = SharedArray{Float64}((nPrimVals,nEvoVals,nSecVals,nRepets))
            SDmeans[category][:] .= NaN
            SDstds[category][:] .= NaN
        end

        loop_vars = [(primInd,rPrimExt,secInd,rSecExt,evoInd,rEvo,repetition) for (primInd,rPrimExt) in enumerate(primVals)
            for (secInd,rSecExt) in enumerate(secVals) for (evoInd,rEvo) in enumerate(evoVals) for repetition in repets];

        @sync @distributed for (primInd,rPrimExt,secInd,rSecExt,evoInd,rEvo,repetition) in loop_vars
            println("$((primInd,rPrimExt,secInd,rSecExt,evoInd,rEvo,repetition))")
            fileName = "data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rSecExt),$(evoInd))_repet=$repetition.jld2" 
            if isfile(fileName)
                sd = load(fileName, "simulationData")
            else
                runFinished[primInd,evoInd,secInd,repetition] = false  #maybe redundant as SharedArray defaults to zeo or here false but not sure if that behaviour is guaranteed
                continue
            end
            if sd.runFinished    #did the run end prematurely?
                SDsPool,specIdsPool = modifiedSDSimilarities(poolnet);
                SDsCol,specIdsCol = modifiedSDSimilarities(sd.colnet);
                SDsColInPool,specIdsColInPool = modifiedSDSimilarities(sd.poolnet,collect(sd.colnet.spec));
                SDsPoolWithoutCol,specIdsPoolWithoutCol = modifiedSDSimilarities(sd.poolnet,collect(setdiff(sd.poolnet.spec,sd.colnet.spec)));
                SDsPoolOnlyEats,specIdsPoolOnlyEats = modifiedSDSimilarities(sd.poolnet,:eat);
                SDsColOnlyEats,specIdsColOnlyEats = modifiedSDSimilarities(sd.colnet,:eat);
                SDsColInPoolOnlyEats,specIdsColInPoolOnlyEats = modifiedSDSimilarities(sd.poolnet,collect(sd.colnet.spec),:eat);
                SDsPoolWithoutColOnlyEats,specIdsPoolWithoutColOnlyEats = modifiedSDSimilarities(sd.poolnet,collect(setdiff(sd.poolnet.spec,sd.colnet.spec)),:eat);
                
                SDmeans[:SDsPool] = nanmean(SDsPool)
                SDmeans[:SDsPoolOnlyEats] = nanmean(SDsPoolOnlyEats)
                SDmeans[:SDsCol] = nanmean(SDsCol)
                SDmeans[:SDsColOnlyEats] = nanmean(SDsColOnlyEats)
                SDmeans[:SDsColInPool] = nanmean(SDsColInPool)
                SDmeans[:SDsColInPoolOnlyEats] = nanmean(SDsColInPoolOnlyEats)
                SDmeans[:SDsPoolWithoutCol] = nanmean(SDsPoolWithoutCol)
                SDmeans[:SDsPoolWithoutColOnlyEats] = nanmean(SDsPoolWithoutColOnlyEats)
                SDstds[:SDsPool] = nanstd(SDsPool)
                SDstds[:SDsPoolOnlyEats] = nanstd(SDsPoolOnlyEats)
                SDstds[:SDsCol] = nanstd(SDsCol)
                SDstds[:SDsColOnlyEats] = nanstd(SDsColOnlyEats)
                SDstds[:SDsColInPool] = nanstd(SDsColInPool)
                SDstds[:SDsColInPoolOnlyEats] = nanstd(SDsColInPoolOnlyEats)
                SDstds[:SDsPoolWithoutCol] = nanstd(SDsPoolWithoutCol)
                SDstds[:SDsPoolWithoutColOnlyEats] = nanstd(SDsPoolWithoutColOnlyEats)
            end
        end
        jldsave("data/$(simulationName)/similarityResults.jld2", compress; SDmeans, SDstds)
    else
        SDmeans, SDstds = 
            load("data/$(simulationName)/similarityResults.jld2", "SDmeans", "SDstds")
    end
    
    plotlyjs()
    xLabel = "evolutionary rate ($(1/couplingFactor)*global extinction rate)"
    seriesLabels = reshape(["rSecExt = $rSecExt" for rSecExt in secVals], (1,nSecVals))
    for (primInd,rPrimExt) in enumerate(primVals)
        plotPath = "data/$simulationName/plots/rPrimExt=$(rPrimExt)"
        mkpath(plotPath)

        for category in categories
            meanPlot = plot(evoVals, dropdims(mean(selectdim(SDmeans[category],1,primInd),dims=3),dims=3),
                xlabel = xLabel, ylabel = "Mean Soerensen-Dice Similarity Coefficient",
                title = "Mean Sorensen-Dice Similarity Coefficient for all $(categoryLongNames[category]).",
                xaxis = :log, labels = seriesLabels, legend = :topright)

            Plots.savefig(meanPlot,"$(plotPath)/meanSD_$(String(category))Plot.svg");

            if any(isnan,SDmeans[category])
                meanPlot = plot(evoVals, NaNStatistics.nanmean(selectdim(SDmeans[category],1,primInd),dim=3),
                xlabel = xLabel, ylabel = "Mean Soerensen-Dice Similarity Coefficient",
                title = "Mean Sorensen-Dice Similarity Coefficient for all $(categoryLongNames[category]) - ignoring NaNs.", 
                xaxis = :log, labels = seriesLabels, legend = :topleft)

                Plots.savefig(maxPlot,"$(plotPath)/meanSD_$(String(vecZParam))Plot_ignoreNaNs.svg");
            end

            stdPlot = plot(evoVals, dropdims(mean(selectdim(SDmeans[category],1,primInd),dims=3),dims=3),
                xlabel = xLabel, ylabel = "Standard deviation of Soerensen-Dice Similarity Coefficient",
                title = "Standard deviation of Sorensen-Dice Similarity Coefficient for all $(categoryLongNames[category]).",
                xaxis = :log, labels = seriesLabels, legend = :topright)

            Plots.savefig(stdPlot,"$(plotPath)/stdSD_$(String(category))Plot.svg");

            if any(isnan,SDstds[category])
                stdPlot = plot(evoVals, NaNStatistics.nanmean(selectdim(SDmeans[category],1,primInd),dim=3),
                    xlabel = xLabel, ylabel = "Standard deviation of Soerensen-Dice Similarity Coefficient",
                    title = "Standard deviation of Sorensen-Dice Similarity Coefficient for all $(categoryLongNames[category]) - ignoring NaNs.",  
                    xaxis = :log, labels = seriesLabels, legend = :topleft)

                Plots.savefig(stdPlot,"$(plotPath)/stdSD_$(String(vecZParam))Plot_ignoreNaNs.svg");
            end
        end
    end
end