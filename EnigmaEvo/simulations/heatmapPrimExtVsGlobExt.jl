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
    paramName = "(rPrimExt,rExt)"
    simulationName = "heatmapPrimExtVsGlobExt";        #specify the name of the simulation
    mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder including plots subfolder
    mkpath("data/$(simulationName)/runs"); 
    compress::Bool = true;  #should the data be compressed before storing?

    jldsave("data/$(simulationName)/parameters.jld2",compress;
        rates0, maxits, cm,ce,cpred, diverse, restrict_colonization,
        logging,S,lambda,SSprobs,SOprobs,nBasalRes)

    primVals = [#=.1,.5,=#1,1.5,2,3,5,10,15,20,40]       #specify the parameter values that shall be simulated
    globVals = [.022,.03,.035,.045,.06,.085,.11] #[.001,.005,.01,.016,.022,.03,.035,.045,.06,.085,.11,.3,.35,.45,.7,1.,2.,5.]
    numPrimVals = length(primVals)
    numGlobVals = length(globVals)

    nRepets = 50
    repets = 1:nRepets

    endVecZParams = [:trophLevels]
    allTimeScalarZParams = Symbol[:specRich, :pool, :meanEats, :meanNeeds,  :nPrimExtSpec, :nSecExtSpec, :nColonizers]
    zParams = [endVecZParams; allTimeScalarZParams]
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
        for vecZParam in endVecZParams
            heatMapsMax[vecZParam] = SharedArray{Float64}((numPrimVals,numGlobVals,nRepets))
            heatMapsMax[vecZParam][:] .= NaN
        end
        for zParam in zParams
            heatMapsMean[zParam] = SharedArray{Float64}((numPrimVals,numGlobVals,nRepets))
            heatMapsMean[zParam][:] .= NaN
        end

        runFinished = SharedArray{Bool}((numPrimVals,numGlobVals,nRepets))

        loop_vars = [(primInd,rPrimExt,globExtInd,rGlobExt,repetition) for (primInd,rPrimExt) in enumerate(primVals)
        for (globExtInd,rGlobExt) in enumerate(globVals) for repetition in repets];

        @sync @distributed for (primInd,rPrimExt,globExtInd,rGlobExt,repetition) in loop_vars
            if onlyPlots
                fileName = "data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rGlobExt))_repet=$repetition.jld2" 
                if isfile(fileName)
                    sd = load(fileName, "simulationData")
                else
                    runFinished[primInd,globExtInd,repetition] = false  #maybe redundant as SharedArray defaults to zeo or here false but not sure if that behaviour is guaranteed
                    continue
                end
            else
                initpoolnet::ENIgMaGraph = setUpPool(S,lambda,nBasalRes,SSprobs,SOprobs,diverse);
                rates0 = (rc = 1., rprimext = rPrimExt, rsecext = 10., reo = 1., revo = 0.35, rext = rGlobExt);
                # EVOLUTIONARY VERSION
                (sd,_) = assemblyevo(initpoolnet, rates0, maxits, cn,cm,ce,cpred, diverse, restrict_colonization);

                jldsave("data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rGlobExt))_repet=$repetition.jld2",compress;
                    simulationData = sd, rates0)
            end
            if 1 <= sd.specRich[end] <= 1000#sd.runFinished    #did the run end prematurely?
                runFinished[primInd,globExtInd,repetition] = true
                for endZParamSymb in endVecZParams
                    endZParam = getfield(sd,endZParamSymb)
                    heatMapsMax[endZParamSymb][primInd,globExtInd,repetition] = mean(maximum.(endZParam))
                    heatMapsMean[endZParamSymb][primInd,globExtInd,repetition] = mean(mean.(endZParam))
                end
                for allTimeZParamSymb in allTimeScalarZParams
                    allTimeZParam = getfield(sd,allTimeZParamSymb)[9500:end]
                    heatMapsMean[allTimeZParamSymb][primInd,globExtInd,repetition] = mean(allTimeZParam)
                end
            else
                runFinished[primInd,globExtInd,repetition] = false
            end
        end

        jldsave("data/$(simulationName)/results.jld2",compress; heatMapsMax, heatMapsMean,runFinished)
    else
        heatMapsMax, heatMapsMean, runFinished = 
            load("data/$(simulationName)/results.jld2", "heatMapsMax", "heatMapsMean", "runFinished")
    end
    
    plotlyjs()

    nFinishedRunsPlot = Plots.heatmap(string.(globVals), string.(primVals),dropdims(sum(runFinished,dims=3),dims=3),
    size = (720,720), xlabel = "global extinction rate", ylabel = "primary extinction rate",
    title = "Number of successfull runs used", xticks = :all, yticks = :all)

    Plots.savefig(nFinishedRunsPlot,"data/$simulationName/plots/nFinishedRunsPlot.html");


    for vecZParam in endVecZParams
        maxPlot = Plots.heatmap(string.(globVals), string.(primVals),dropdims(mean(heatMapsMax[vecZParam],dims=3),dims=3),
        size = (720,720), xlabel = "global extinction rate", ylabel = "primary extinction rate",
        title = "Average maximal $(zParamLongNames[vecZParam]) in itterations 9500 to 10000", xticks = :all, yticks = :all)

        Plots.savefig(maxPlot,"data/$simulationName/plots/max_$(String(vecZParam))Plot.html");

        if any(isnan,heatMapsMax[vecZParam])
            maxPlot = Plots.heatmap(string.(globVals), string.(primVals),NaNStatistics.nanmean(heatMapsMax[vecZParam],dim=3),
            size = (720,720), xlabel = "global extinction rate", ylabel = "primary extinction rate",
            title = "Average maximal $(zParamLongNames[vecZParam]) in itterations 9500 to 10000 - ignoring NaNs.", xticks = :all, yticks = :all)

            Plots.savefig(maxPlot,"data/$simulationName/plots/max_$(String(vecZParam))Plot_ignoreNaNs.html");
        end
    end

    for zParam in zParams
        meanPlot = Plots.heatmap(string.(globVals), string.(primVals),dropdims(mean(heatMapsMean[zParam],dims=3),dims=3),
            size = (720,720), xlabel = "global extinction rate", ylabel = "primary extinction rate",
            title = "Average mean $(zParamLongNames[zParam]) in itterations 9500 to 10000", xticks = :all, yticks = :all)

        Plots.savefig(meanPlot,"data/$simulationName/plots/mean_$(String(zParam))Plot.html");

        if any(isnan,heatMapsMean[zParam])
            meanPlot = Plots.heatmap(string.(globVals), string.(primVals),NaNStatistics.nanmean(heatMapsMean[zParam],dim=3),
            size = (720,720), xlabel = "global extinction rate", ylabel = "primary extinction rate",
            title = "Average mean $(zParamLongNames[zParam]) in itterations 9500 to 10000 - ignoring NaNs", xticks = :all, yticks = :all)

            Plots.savefig(meanPlot,"data/$simulationName/plots/mean_$(String(zParam))Plot_ignoreNaNs.html");
        end
    end
end

#simulation()

#=function createPlotsFromFiles(;fromResultsFile = false)   #supposedly its better to wrap stuff in functions and not use global variables if possible(not sure if typed ones are a problem)
    include("src/set_up_params.jl");
    maxits = 10_000
    #prepare everything for a simulation consisting of the variation of a parmeter
    paramName = "(rPrimExt,rExt)"
    simulationName = "heatmapPrimExtVsGlobExt";        #specify the name of the simulation
    mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder including plots subfolder
    mkpath("data/$(simulationName)/runs"); 
    compress::Bool = true;  #should the data be compressed before storing?

    jldsave("data/$(simulationName)/parameters.jld2",compress;
        rates0, maxits, cm,ce,cpred, diverse, restrict_colonization,
        logging,S,lambda,SSprobs,SOprobs,nBasalRes)

    primVals = [.1,.5,1,1.5,2,3,5,10,15,20,40]       #specify the parameter values that shall be simulated
    globVals = [.001,.005,.01,.016,.3,.35,.45,.7,1.,2.,5.]
    numPrimVals = length(primVals)
    numGlobVals = length(globVals)

    nRepets = 50
    repets = 1:nRepets

    endZParams = [:trophLevels]
    allTimeZParams = Symbol[:specRich, :pool, :meanEats, :meanNeeds,  :nPrimExtSpec, :nSecExtSpec, :nColonizers]
    zParams = [endZParams; allTimeZParams]
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

    if !fromResultsFile
    heatMapsMax = Dict{Symbol, SharedArray{Float64}}()
    heatMapsMean = Dict{Symbol,SharedArray{Float64}}()
    for zParam in zParams
        heatMapsMax[zParam] = SharedArray{Float64}((numPrimVals,numGlobVals,nRepets))
        heatMapsMean[zParam] = SharedArray{Float64}((numPrimVals,numGlobVals,nRepets))
        heatMapsMax[zParam][:] .= NaN
        heatMapsMean[zParam][:] .= NaN
    end

    runFinished = SharedArray{Bool}((numPrimVals,numGlobVals,nRepets))

    loop_vars = [(primInd,rPrimExt,globExtInd,rGlobExt,repetition) for (primInd,rPrimExt) in enumerate(primVals)
    for (globExtInd,rGlobExt) in enumerate(globVals) for repetition in repets];

    @sync @distributed for (primInd,rPrimExt,globExtInd,rGlobExt,repetition) in loop_vars
        fileName = "data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rGlobExt))_repet=$repetition.jld2" 
        if isfile(fileName)
            sd = load(fileName, "simulationData")
            if 1 <= sd.specRich[end] <= maxits    #did the run end prematurely (no valid values at end) not perfect way to check that
                runFinished[primInd,globExtInd,repetition] = true
                for endZParamSymb in endZParams
                    endZParam = getfield(sd,endZParamSymb)
                    heatMapsMax[endZParamSymb][primInd,globExtInd,repetition] = mean(maximum.(endZParam))
                    heatMapsMean[endZParamSymb][primInd,globExtInd,repetition] = mean(mean.(endZParam))
                end
                for allTimeZParamSymb in allTimeZParams
                    allTimeZParam = getfield(sd,allTimeZParamSymb)[9500:end]
                    heatMapsMax[allTimeZParamSymb][primInd,globExtInd,repetition] = mean(maximum.(allTimeZParam))
                    heatMapsMean[allTimeZParamSymb][primInd,globExtInd,repetition] = mean(mean.(allTimeZParam))
                end
            else
                runFinished[primInd,globExtInd,repetition] = false
            end
        end
    end

        jldsave("data/$(simulationName)/results.jld2",compress; heatMapsMax, heatMapsMean,runFinished)
    else
        heatMapsMax, heatMapsMean = load("data/$(simulationName)/results.jld2", "heatMapsMax", "heatMapsMean")
    end
    
    plotlyjs()

    #nFinishedRunsPlot = Plots.heatmap(primVals, globVals, dropdims(sum(runFinished,dims=3),dims=3),
    #size = (1280,720), xlabel = "primary extinction rate", ylabel = "global extinction rate",
    #title = "Number of successfull runs used")


    #Plots.savefig(nFinishedRunsPlot,"data/$simulationName/plots/nFinishedRunsPlot.html");

    @distributed for zParam in zParams
        maxPlot = Plots.heatmap(primVals, globVals, dropdims(mean(heatMapsMax[zParam],dims=3),dims=3),
            size = (1280,720), xlabel = "global extinction rate", ylabel = "primary extinction rate",
            title = "Average maximal $(zParamLongNames[zParam]) in itterations 9500 to 10000")
        meanPlot = Plots.heatmap(primVals, globVals, dropdims(mean(heatMapsMean[zParam],dims=3),dims=3),
            size = (1280,720), xlabel = "primary extinction rate", ylabel = "global extinction rate",
            title = "Average mean $(zParamLongNames[zParam]) in itterations 9500 to 10000")


        Plots.savefig(maxPlot,"data/$simulationName/plots/max_$(String(zParam))Plot.html");
        Plots.savefig(meanPlot,"data/$simulationName/plots/mean_$(String(zParam))Plot.html");

        if any(isnan,heatMapsMax[zParam])
            maxPlot = Plots.heatmap(primVals, globVals, NaNStatistics.nanmean(heatMapsMax[zParam],dim=3),
            size = (1280,720), xlabel = "global extinction rate", ylabel = "primary extinction rate",
            title = "Average maximal $(zParamLongNames[zParam]) in itterations 9500 to 10000 - ignoring NaNs.")

            Plots.savefig(maxPlot,"data/$simulationName/plots/max_$(String(zParam))Plot_ignoreNaNs.html");
        end

        if any(isnan, heatMapsMean[zParam])
            meanPlot = Plots.heatmap(primVals, globVals, NaNStatistics.nanmean(heatMapsMean[zParam],dim=3),
            size = (1280,720), xlabel = "primary extinction rate", ylabel = "global extinction rate",
            title = "Average mean $(zParamLongNames[zParam]) in itterations 9500 to 10000 - ignoring NaNs")

            Plots.savefig(meanPlot,"data/$simulationName/plots/mean_$(String(zParam))Plot_ignoreNaNs.html");
        end
    end
end=#