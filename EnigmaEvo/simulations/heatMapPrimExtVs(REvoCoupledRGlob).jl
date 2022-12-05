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

function simulation_primExtVsREvoCoupledRGlob(onlyPlots=false,fromResults=false;couplingFactor=.5)   #supposedly its better to wrap stuff in functions and not use global variables if possible(not sure if typed ones are a problem)
    if fromResults && !onlyPlots
        error("If only results file should be used it implies that only plots can be created.")
    end
    include("src/set_up_params.jl");
    maxits = 10_000
    #prepare everything for a simulation consisting of the variation of a parmeter
    paramName = "(rPrimExt,rEvo)"
    simulationName = "heatmapPrimExtVs(REvoCoupledRGlob)";        #specify the name of the simulation
    mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder including plots subfolder
    mkpath("data/$(simulationName)/runs"); 
    compress::Bool = true;  #should the data be compressed before storing?
    
    primVals = [1,1.5,2,2.5,3,4,5,10,15,20]       #specify the parameter values that shall be simulated
    rEvoVals = [0.005, 0.01, 0.016, 0.022, 0.03, 0.035, 0.045, 0.06, 0.08, 0.22, 0.3, 0.35, 0.45, 0.55, 0.7, 1.0, 2.0]#0.015:.0005:0.022
    numPrimVals = length(primVals)
    numREvoVals = length(rEvoVals)
    nRepets = 50
    repets = 1:nRepets

    jldsave("data/$(simulationName)/parameters.jld2",compress;
        rates0, maxits, cm,ce,cpred, diverse, restrict_colonization,
        logging,S,lambda,SSprobs,SOprobs,nBasalRes,primVals, rEvoVals, nRepets)
        
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
            heatMapsMax[vecZParam] = SharedArray{Float64}((numPrimVals,numREvoVals,nRepets))
            heatMapsMax[vecZParam][:] .= NaN
        end
        for zParam in zParams
            heatMapsMean[zParam] = SharedArray{Float64}((numPrimVals,numREvoVals,nRepets))
            heatMapsMean[zParam][:] .= NaN
        end

        runFinished = SharedArray{Bool}((numPrimVals,numREvoVals,nRepets))

        loop_vars = [(primInd,rPrimExt,rEvoInd,rEvo,repetition) for repetition in repets
            for (primInd,rPrimExt) in enumerate(primVals) for (rEvoInd,rEvo) in enumerate(rEvoVals)];

        @sync @distributed for (primInd,rPrimExt,rEvoInd,rEvo,repetition) in loop_vars
            if onlyPlots
                fileName = "data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rEvo))_repet=$repetition.jld2" 
                if isfile(fileName)
                    sd = load(fileName, "simulationData")
                else
                    runFinished[primInd,rEvoInd,repetition] = false  #maybe redundant as SharedArray defaults to zeo or here false but not sure if that behaviour is guaranteed
                    continue
                end
            else
                initpoolnet::ENIgMaGraph = setUpPool(S,lambda,nBasalRes,SSprobs,SOprobs,diverse);
                rates0 = (rc = 1., rprimext = rPrimExt, rsecext = 10., reo = 1., revo = rEvo, rext = couplingFactor*rEvo);
                # EVOLUTIONARY VERSION
                (sd,_) = assemblyevo(initpoolnet, rates0, maxits, cn,cm,ce,cpred, diverse, restrict_colonization);

                jldsave("data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rEvo))_repet=$repetition.jld2",compress;
                    simulationData = sd, rates0)
            end
            if sd.runFinished    #did the run end prematurely?
                runFinished[primInd,rEvoInd,repetition] = true
                for endZParamSymb in endVecZParams
                    endZParam = getfield(sd,endZParamSymb)
                    heatMapsMax[endZParamSymb][primInd,rEvoInd,repetition] = mean(maximum.(endZParam))
                    heatMapsMean[endZParamSymb][primInd,rEvoInd,repetition] = mean(mean.(endZParam))
                end
                for allTimeZParamSymb in allTimeScalarZParams
                    allTimeZParam = getfield(sd,allTimeZParamSymb)[9500:end]
                    heatMapsMean[allTimeZParamSymb][primInd,rEvoInd,repetition] = mean(allTimeZParam)
                end
            else
                runFinished[primInd,rEvoInd,repetition] = false
            end
        end

        jldsave("data/$(simulationName)/results.jld2",compress; heatMapsMax, heatMapsMean,runFinished)
    else
        heatMapsMax, heatMapsMean, runFinished = 
            load("data/$(simulationName)/results.jld2", "heatMapsMax", "heatMapsMean", "runFinished")
    end
    
    plotlyjs()

    nFinishedRunsPlot = Plots.heatmap("(".*string.(rEvoVals).*", ".*string.(rEvoVals.*couplingFactor).*")", string.(primVals), dropdims(sum(runFinished,dims=3),dims=3),
        size = (720,720), xlabel = "(evolutionary rate, global extinction rate)", ylabel = "primary extinction rate",
        title = "Number of successfull runs used", xrotation = 60, xticks = :all, yticks = :all)

    Plots.savefig(nFinishedRunsPlot,"data/$simulationName/plots/nFinishedRunsPlot.html");


    for vecZParam in endVecZParams
        maxPlot = Plots.heatmap("(".*string.(rEvoVals).*", ".*string.(rEvoVals.*couplingFactor).*")", string.(primVals), dropdims(mean(heatMapsMax[vecZParam],dims=3),dims=3),
            size = (720,720), xlabel = "(evolutionary rate, global extinction rate)", ylabel = "primary extinction rate",
            title = "Average maximal $(zParamLongNames[vecZParam]) in itterations 9500 to 10000", xrotation = 60, xticks = :all, yticks = :all)

        Plots.savefig(maxPlot,"data/$simulationName/plots/max_$(String(vecZParam))Plot.html");

        if any(isnan,heatMapsMax[vecZParam])
            maxPlot = Plots.heatmap("(".*string.(rEvoVals).*", ".*string.(rEvoVals.*couplingFactor).*")", string.(primVals), NaNStatistics.nanmean(heatMapsMax[vecZParam],dim=3),
                size = (720,720), xlabel = "(evolutionary rate, global extinction rate)", ylabel = "primary extinction rate",
                title = "Average maximal $(zParamLongNames[vecZParam]) in itterations 9500 to 10000 - ignoring NaNs.", xrotation = 60, xticks = :all, yticks = :all)

            Plots.savefig(maxPlot,"data/$simulationName/plots/max_$(String(vecZParam))Plot_ignoreNaNs.html");
        end
    end

    for zParam in zParams
        meanPlot = Plots.heatmap("(".*string.(rEvoVals).*", ".*string.(rEvoVals.*couplingFactor).*")", string.(primVals), dropdims(mean(heatMapsMean[zParam],dims=3),dims=3),
            size = (720,720), xlabel = "(evolutionary rate, global extinction rate)", ylabel = "primary extinction rate",
            title = "Average mean $(zParamLongNames[zParam]) in itterations 9500 to 10000", xrotation = 60, xticks = :all, yticks = :all)

        Plots.savefig(meanPlot,"data/$simulationName/plots/mean_$(String(zParam))Plot.html");

        if any(isnan,heatMapsMean[zParam])
            meanPlot = Plots.heatmap("(".*string.(rEvoVals).*", ".*string.(rEvoVals.*couplingFactor).*")", string.(primVals), NaNStatistics.nanmean(heatMapsMean[zParam],dim=3),
                size = (720,720), xlabel = "(evolutionary rate, global extinction rate)", ylabel = "primary extinction rate",
                title = "Average mean $(zParamLongNames[zParam]) in itterations 9500 to 10000 - ignoring NaNs", xrotation = 60, xticks = :all, yticks = :all)

            Plots.savefig(meanPlot,"data/$simulationName/plots/mean_$(String(zParam))Plot_ignoreNaNs.html");
        end
    end
end     