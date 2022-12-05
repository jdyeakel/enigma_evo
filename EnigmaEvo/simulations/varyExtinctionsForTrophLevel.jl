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

function simulation(onlyPlots=false,fromResults=false)   #supposedly its better to wrap stuff in functions and not use global variables if possible(not sure if typed ones are a problem)
    if fromResults && !onlyPlots
        error("If only results file should be used it implies that only plots can be created.")
    end

    include("src/set_up_params.jl");
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

    jldsave("data/$(simulationName)/parameters.jld2",compress;
    rates0, maxits, cm,ce,cpred, diverse, restrict_colonization,
    logging,S,lambda,SSprobs,SOprobs,nBasalRes,paramVals,nRepets)

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
        heatMapsMax = Dict{Symbol, SharedArray{Float64}}()
        heatMapsMean = Dict{Symbol,SharedArray{Float64}}()
        for vecZParam in endVectorialZParams
            heatMapsMax[vecZParam] = SharedArray{Float64}((numParams,numParams,nRepets))
            heatMapsMax[vecZParam][:] .= NaN
        end
        for zParam in zParams
            heatMapsMean[zParam] = SharedArray{Float64}((numParams,numParams,nRepets))
            heatMapsMean[zParam][:] .= NaN
        end

        runFinished = SharedArray{Bool}((numParams,numParams,nRepets))

        loop_vars = [(primInd,rPrimExt,secInd,rSecExt,repetition) for (primInd,rPrimExt) in enumerate(paramVals)
        for (secInd,rSecExt) in enumerate(paramVals) for repetition in repets];

        @sync @distributed for (primInd,rPrimExt,secInd,rSecExt,repetition) in loop_vars
            if onlyPlots
                fileName = "data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rSecExt))_repet=$repetition.jld2" 
                if isfile(fileName)
                    sd = load(fileName, "simulationData")
                else
                    runFinished[primInd,secInd,repetition] = false  #maybe redundant as SharedArray defaults to zeo or here false but not sure if that behaviour is guaranteed
                    continue
                end
            else
                initpoolnet::ENIgMaGraph = setUpPool(S,lambda,nBasalRes,SSprobs,SOprobs,diverse);
                rates0 = (rc = 1., rprimext = rPrimExt, rsecext = rSecExt, reo = 1., revo = 0.35, rext = 0.16);
                # EVOLUTIONARY VERSION
                (sd,_) = assemblyevo(initpoolnet, rates0, maxits, cn,cm,ce,cpred, diverse, restrict_colonization);

                jldsave("data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rSecExt))_repet=$repetition.jld2",compress;
                    simulationData = sd, rates0)
            end
            if sd.runFinished    #did the run end prematurely?
                runFinished[primInd,secInd,repetition] = true
                for endvectorialZParamSymb in endVectorialZParams
                    endVectorialZParam = getfield(sd,endvectorialZParamSymb)
                    heatMapsMax[endvectorialZParamSymb][primInd,secInd,repetition] = mean(maximum.(endVectorialZParam))
                    heatMapsMean[endvectorialZParamSymb][primInd,secInd,repetition] = mean(mean.(endVectorialZParam))
                end
                for allTimeScalarZParamSymb in allTimeScalarZParams
                    allTimeScalarZParam = getfield(sd,allTimeScalarZParamSymb)[9500:end]
                    heatMapsMean[allTimeScalarZParamSymb][primInd,secInd,repetition] = mean(allTimeScalarZParam)
                end
            else
                runFinished[primInd,secInd,repetition] = false
            end
        end

        jldsave("data/$(simulationName)/results.jld2",compress; heatMapsMax, heatMapsMean,runFinished)
    else
        heatMapsMax, heatMapsMean, runFinished = 
        load("data/$(simulationName)/results.jld2", "heatMapsMax", "heatMapsMean", "runFinished")
    end
    
    plotlyjs()

    nFinishedRunsPlot = Plots.heatmap(string.(paramVals), string.(paramVals),  dropdims(sum(runFinished,dims=3),dims=3),
    size = (720,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
    title = "Number of successfull runs used", xticks = :all, yticks = :all, xrotation = 60)

    Plots.savefig(nFinishedRunsPlot,"data/$simulationName/plots/nFinishedRunsPlot.html");

    for vecZParam in endVectorialZParams
        maxPlot = Plots.heatmap(string.(paramVals), string.(paramVals), dropdims(mean(heatMapsMax[vecZParam],dims=3),dims=3),
        size = (720,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
        title = "Average maximal $(zParamLongNames[vecZParam]) in itterations 9500 to 10000", xticks = :all,
        yticks = :all, xrotation = 60)

        Plots.savefig(maxPlot,"data/$simulationName/plots/max_$(String(vecZParam))Plot.html");

        if any(isnan,heatMapsMax[vecZParam])
            maxPlot = Plots.heatmap(string.(paramVals), string.(paramVals), NaNStatistics.nanmean(heatMapsMax[vecZParam],dim=3),
            size = (720,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
            title = "Average maximal $(zParamLongNames[vecZParam]) in itterations 9500 to 10000 - ignoring NaNs.", 
            xticks = :all, yticks = :all, xrotation = 60)

            Plots.savefig(maxPlot,"data/$simulationName/plots/max_$(String(vecZParam))Plot_ignoreNaNs.html");
        end
    end

    for zParam in zParams
        meanPlot = Plots.heatmap(string.(paramVals), string.(paramVals), dropdims(mean(heatMapsMean[zParam],dims=3),dims=3),
            size = (720,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
            title = "Average mean $(zParamLongNames[zParam]) in itterations 9500 to 10000", xticks = :all,
            yticks = :all, xrotation = 60)

        Plots.savefig(meanPlot,"data/$simulationName/plots/mean_$(String(zParam))Plot.html");

        if any(isnan,heatMapsMean[zParam])
            meanPlot = Plots.heatmap(string.(paramVals), string.(paramVals), NaNStatistics.nanmean(heatMapsMean[zParam],dim=3),
            size = (720,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
            title = "Average mean $(zParamLongNames[zParam]) in itterations 9500 to 10000 - ignoring NaNs",
            xticks = :all, yticks = :all, xrotation = 60)

            Plots.savefig(meanPlot,"data/$simulationName/plots/mean_$(String(zParam))Plot_ignoreNaNs.html");
        end
    end
end