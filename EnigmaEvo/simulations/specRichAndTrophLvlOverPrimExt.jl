if !@isdefined(NO_STANDALONE) || NO_STANDALONE == false
    if homedir() == "/home/z840"    #for downward compatibility ;)
        localPath::String = "$(homedir())/2019_Lego_Evo/EnigmaEvo/src/";
    elseif isfile("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl")
        localPath::String = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/";
    else                            #othserwise use relative paths
        localPath::String = "../src/";
    end;

    #using Distributed
    #addprocs(80)

    include(localPath*"loadfuncs.jl");

    using SharedArrays
    @everywhere import NaNStatistics
end

""""
    rPrimExtLinePlot(onlyPlots=false, fromResults=false; plotLib = :Plots)

    Performs a series of simulations varying over the primary extinction rate. Creates line plots of results.

    # Arguments
    'onlyPlots::Bool': If true assumes the simulations to be finished and only performs data analysis from files and plots. 
    'fromResults::Bool': If true loads final results and plots. Overrides 'onlyPlots'.
    'plotLib::Symbol': Symbolic name of plotting library. Supported are :Plots and :Makie
"""
function rPrimExtLinePlot(onlyPlots=false,fromResults=false; plotLib = :Plots)   #supposedly its better to wrap stuff in functions and not use global variables if possible(not sure if typed ones are a problem)
    if fromResults && !onlyPlots
        error("If only results file should be used it implies that only plots can be created.")
    end

    include("src/set_up_params.jl");
    maxits = 30_000
    #prepare everything for a simulation consisting of the variation of a parmeter
    paramName = "rPrimExt"
    simulationName = "specRichAndTrophLvlOverPrimExt2";        #specify the name of the simulation
    mkpath("data/$(simulationName)/runs");    #make a folder with that name in the Data folder including plots subfolder
    mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder including plots subfolder
    compress::Bool = true;  #should the data be compressed before storing?


    primVals = 1.:.5:20.       #specify the parameter values that shall be simulated
    nPrimVals = length(primVals)
    secVals = [1.,10.]
    nSecVals = length(secVals)
    nRepets = 50
    repets = 1:nRepets

    jldsave("data/$(simulationName)/parameters.jld2",compress;
    rates0, maxits, cm,ce,cpred, diverse, restrict_colonization,
    logging,S,lambda,SSprobs,SOprobs,nBasalRes,primVals,nRepets, secVals)

    #define names of properties that are analysed. Somewhat messy categorized by the way they are collected
    endVectorialZParams = [:trophLevels] # properties that are only computed at the last simulations (to save time) and are vectorial (here value for each species)
    allTimeScalarZParams = Symbol[:specRich, :pool, :meanEats, :meanNeeds,  :nPrimExtSpec, :nSecExtSpec, :nColonizers]  #scalar properties of the whole network collected for every itteration of the run
    zParams = [endVectorialZParams; allTimeScalarZParams]
    zParamLongNames = Dict{Symbol,String}(
        :trophLevels => "trophic level",
        :specRich => "species richness",
        :pool => "species richness in the pool",
        :meanEats => "mean number of eat interactions",
        :meanNeeds => "mean number of need interactions",
        :nPrimExtSpec => "number of primary extinction candidates",
        :nSecExtSpec => "number of secondary extinction candidates",
        :nColonizers => "number of potential Colonizers"
    )

    if !fromResults
        maxResults = Dict{Symbol, SharedArray{Float64}}()
        meanResults = Dict{Symbol,SharedArray{Float64}}()
        for vecZParam in endVectorialZParams
            maxResults[vecZParam] = SharedArray{Float64}((nPrimVals,nSecVals,nRepets))
            maxResults[vecZParam][:] .= NaN
        end
        for zParam in zParams
            meanResults[zParam] = SharedArray{Float64}((nPrimVals,nSecVals,nRepets))
            meanResults[zParam][:] .= NaN
        end

        runFinished = SharedArray{Bool}((nPrimVals,nSecVals,nRepets))

        loop_vars = [(primInd,rPrimExt,secInd,rSecExt,repetition) for (primInd,rPrimExt) in enumerate(primVals)
            for (secInd,rSecExt) in enumerate(secVals) for repetition in repets];

        @sync @distributed for (primInd,rPrimExt,secInd,rSecExt,repetition) in loop_vars
            println("$((primInd,rPrimExt,secInd,rSecExt,repetition))")
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
                (sd,_) = assemblyevo(initpoolnet, rates0, maxits, cm,cn,ce,cpred, diverse, restrict_colonization);

                jldsave("data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rSecExt))_repet=$repetition.jld2",compress;
                    simulationData = sd, rates0)
            end
            if sd.runFinished    #did the run end prematurely?
                runFinished[primInd,secInd,repetition] = true
                for endvectorialZParamSymb in endVectorialZParams
                    endVectorialZParam = getfield(sd,endvectorialZParamSymb)
                    if any(isempty,endVectorialZParam)
                        maxResults[endvectorialZParamSymb][primInd,secInd,repetition] = NaN
                        meanResults[endvectorialZParamSymb][primInd,secInd,repetition] = NaN
                    else
                        maxResults[endvectorialZParamSymb][primInd,secInd,repetition] = mean(maximum.(endVectorialZParam))
                        meanResults[endvectorialZParamSymb][primInd,secInd,repetition] = mean(mean.(endVectorialZParam))
                    end
                end
                for allTimeScalarZParamSymb in allTimeScalarZParams
                    allTimeScalarZParam = getfield(sd,allTimeScalarZParamSymb)[(maxits-1500):end]
                    meanResults[allTimeScalarZParamSymb][primInd,secInd,repetition] = mean(allTimeScalarZParam)
                end
            else
                runFinished[primInd,secInd,repetition] = false
            end
        end
        jldsave("data/$(simulationName)/results.jld2",compress; maxResults, meanResults,runFinished)
    else
        maxResults, meanResults, runFinished = 
        load("data/$(simulationName)/results.jld2", "maxResults", "meanResults", "runFinished")
    end
    
    seriesLabels = reshape([L"r_{2\degree} = %$rSecExt" for rSecExt in secVals], (1,nSecVals))

    #select plot library and do plots
    if plotLib == :Plots
        gr()#plotlyjs()

        nFinishedRunsPlot = scatter(primVals,  dropdims(sum(runFinished,dims=3),dims=3),
        xlabel = "primary extinction rate", ylabel = "number of successfull runs used",
        title = "Number of successfull runs used", labels = seriesLabels)

        Plots.savefig(nFinishedRunsPlot,"data/$simulationName/plots/nFinishedRunsPlot.html");

        for vecZParam in endVectorialZParams
            maxPlot = plot(primVals, dropdims(mean(maxResults[vecZParam],dims=3),dims=3),
            xlabel = "primary extinction rate", ylabel = zParamLongNames[vecZParam],
            #title = "Average maximal $(zParamLongNames[vecZParam]) in itterations $(maxits - 1500) to $maxits",
            labels = seriesLabels)

            Plots.savefig(maxPlot,"data/$simulationName/plots/max_$(String(vecZParam))Plot.html");

            if any(isnan,maxResults[vecZParam])
                maxPlot = plot(primVals, NaNStatistics.nanmean(maxResults[vecZParam],dim=3),
                xlabel = "primary extinction rate", ylabel = zParamLongNames[vecZParam],
                #title = "Average maximal $(zParamLongNames[vecZParam]) in itterations $(maxits - 1500) to $maxits - ignoring NaNs.", 
                labels = seriesLabels)

                Plots.savefig(maxPlot,"data/$simulationName/plots/max_$(String(vecZParam))Plot_ignoreNaNs.html");
            end
        end

        for zParam in zParams
            meanPlot = plot(primVals, dropdims(mean(meanResults[zParam],dims=3),dims=3),
                xlabel = "primary extinction rate", ylabel = zParamLongNames[zParam],
                #title = "Average mean $(zParamLongNames[zParam]) in itterations $(maxits - 1500) to $maxits",
                yaxis = :log, labels = seriesLabels)

            Plots.savefig(meanPlot,"data/$simulationName/plots/mean_$(String(zParam))Plot.html");

            if any(isnan,meanResults[zParam])
                meanPlot = plot(primVals, NaNStatistics.nanmean(meanResults[zParam],dim=3),
                xlabel = "primary extinction rate", ylabel = zParamLongNames[zParam],
                #title = "Average mean $(zParamLongNames[zParam]) in itterations $(maxits - 1500) to $maxits - ignoring NaNs",
                yaxis = :log, labels = seriesLabels)

                Plots.savefig(meanPlot,"data/$simulationName/plots/mean_$(String(zParam))Plot_ignoreNaNs.html");
            end
        end

        zParam = :meanEats

        meanPlot = plot(primVals, diff(dropdims(mean(meanResults[zParam],dims=3),dims=3),dims=2),
        xlabel = "primary extinction rate", ylabel = " mean eats for rSecExt = 10 - mean eats for rSecExt = 1")

        Plots.savefig(meanPlot,"data/$simulationName/plots/meanEatsDifferencePlot.html");

        meanEats = dropdims(mean(meanResults[zParam],dims=3),dims=3)
        meanPlot = plot(primVals, diff(meanEats,dims=2),
        xlabel = "primary extinction rate", ylabel = " mean eats for rSecExt = 10 - mean eats for rSecExt = 1")

        Plots.savefig(meanPlot,"data/$simulationName/plots/meanEatsDifferencePlot.html");
    
    elseif plotLib == :Makie
        unboundedMin = 1
        unboundedMax = 3
        axisLimits = ((1,nothing),nothing)
        xLabel = L"r_{1\degree}"
        for zParam in endVectorialZParams
            fig = CM.Figure(resolution = plotSizePt, fontsize = fontSize, figure_padding = figurePadding);
            ax = CM.Axis(fig[1,1];xlabel = xLabel, ylabel = zParamLongNames[zParam],
                        yscale = identity, limits = axisLimits);
            se = CM.series!(ax,primVals,transpose(dropdims(mean(maxResults[zParam],dims=3),dims=3)),
                labels = seriesLabels, linewidth = lineWidth, color = plotColors);
            CM.vspan!(unboundedMin,unboundedMax,color=unboundedColor,label="species richness\nunbounded",labelfont = :bold,overdraw = true)
            CM.axislegend(ax,patchsize=legendPatchSize, position = legendPosition)   
    
            save("data/$simulationName/plots/max_$(String(zParam))Plot.$(plotDataType)", fig, pt_per_unit = 1);

            if any(isnan,maxResults[zParam])
                fig = CM.Figure(resolution = plotSizePt, fontsize = fontSize, figure_padding = figurePadding);
                ax = CM.Axis(fig[1,1];xlabel = xLabel, ylabel = zParamLongNames[zParam],
                            yscale = log10, limits = axisLimits);
                se = CM.series!(ax,primVals,transpose(NaNStatistics.nanmean(maxResults[zParam],dim=3)),
                    labels = seriesLabels, linewidth = lineWidth, color = plotColors);
                CM.vspan!(unboundedMin,unboundedMax,color=unboundedColor,label="species richness\nunbounded",labelfont = :bold,overdraw = true)
                CM.axislegend(ax,patchsize=legendPatchSize, position = legendPosition)   
        
                save("data/$simulationName/plots/max_$(String(zParam))Plot_ignoreNaNs.$(plotDataType)", fig, pt_per_unit = 1);
            end
        end

        for zParam in zParams
            fig = CM.Figure(resolution = plotSizePt, fontsize = fontSize, figure_padding = figurePadding);
            ax = CM.Axis(fig[1,1];xlabel = xLabel, ylabel = zParamLongNames[zParam],
                        yscale = log10, limits = axisLimits);
            se = CM.series!(ax,primVals,transpose(dropdims(mean(meanResults[zParam],dims=3),dims=3)),
                labels = seriesLabels, linewidth = lineWidth, color = plotColors);
            CM.vspan!(unboundedMin,unboundedMax,color=unboundedColor,label="species richness\nunbounded",labelfont = :bold,overdraw = true)
            CM.axislegend(ax,patchsize=legendPatchSize, position = legendPosition)   
    
            save("data/$simulationName/plots/mean_$(String(zParam))Plot.$(plotDataType)", fig, pt_per_unit = 1);

            if any(isnan,meanResults[zParam])
                fig = CM.Figure(resolution = plotSizePt, fontsize = fontSize, figure_padding = figurePadding);
                ax = CM.Axis(fig[1,1];xlabel = xLabel, ylabel = zParamLongNames[zParam],
                            yscale = log10, limits = axisLimits);
                se = CM.series!(ax,primVals,transpose(NaNStatistics.nanmean(meanResults[zParam],dim=3)),
                    labels = seriesLabels, linewidth = lineWidth, color = plotColors);
                CM.vspan!(unboundedMin,unboundedMax,color=unboundedColor,label="species richness\nunbounded",labelfont = :bold,overdraw = true)
                CM.axislegend(ax,patchsize=legendPatchSize, position = legendPosition)   
        
                save("data/$simulationName/plots/mean_$(String(zParam))Plot_ignoreNaNs.$(plotDataType)", fig, pt_per_unit = 1);
            end
        end

        #final plots
        fig = CM.Figure(resolution = doublePlotSizePt, fontsize = fontSize, 
            figure_padding = figurePadding);
        gridLayout = fig.layout
        zParam = :specRich
        axt = CM.Axis(gridLayout[1,1];ylabel = zParamLongNames[zParam],
            yscale = log10, limits = axisLimits);
        CM.hidexdecorations!(axt,grid = false)
        se = CM.series!(primVals,transpose(dropdims(mean(meanResults[zParam],dims=3),dims=3)),
            labels = seriesLabels, linewidth = 1, color = plotColors);
        CM.vspan!(1,3,color=unboundedColor,label="species richness\nunbounded",
            labelfont = :bold,overdraw = true)
            CM.axislegend(axt,patchsize=legendPatchSize, position = legendPosition)   
    
        zParam = :trophLevels
        axb = CM.Axis(gridLayout[2,1];xlabel = L"r_{1\degree}",ylabel = "maximal $(zParamLongNames[zParam])",
            yscale = identity, limits = axisLimits, ylabelpadding = 17.5);
        se = CM.series!(primVals,transpose(dropdims(mean(maxResults[zParam],dims=3),dims=3)),
            labels = seriesLabels, linewidth = 1, color = plotColors);
        CM.vspan!(1,3,color=unboundedColor,
            label="species richness\nunbounded",labelfont = :bold,overdraw = true)
        
        CM.rowgap!(gridLayout,10)
        CM.linkxaxes!(axb,axt)
    
        for (label, layout) in zip(["A", "B"], [fig[1,1], fig[2,1]])
            CM.Label(layout[1, 1, CM.TopLeft()], label,
                #fontsize = 26,
                font = :bold,
                padding = (0, 29, 5, 0),
                halign = :left)
        end

        save("data/finalPlots/primExtSpecTroph.$(plotDataType)", fig, pt_per_unit = 1);

        fig = CM.Figure(resolution = doublePlotSizePt, fontsize = fontSize, figure_padding = figurePadding);
        gridLayout = fig.layout
        zParam = :meanEats
        axt = CM.Axis(gridLayout[1,1];ylabel = zParamLongNames[zParam],
                            yscale = identity, limits = axisLimits, ylabelpadding = 7.3);
        CM.hidexdecorations!(axt,grid = false)
        se = CM.series!(primVals,transpose(dropdims(mean(meanResults[zParam],dims=3),dims=3)),
            labels = seriesLabels, linewidth = 1, color = plotColors);
        CM.vspan!(1,3,color=unboundedColor,label="species richness\nunbounded",
            labelfont = :bold,overdraw = true)
        CM.axislegend(axt,patchsize=(20,10))
    
        zParam = :meanNeeds
        axb = CM.Axis(gridLayout[2,1]; xlabel = L"r_{1\degree}", ylabel = zParamLongNames[zParam],
            yscale = identity, limits = axisLimits)#, ylabelpadding = 17.5);
        se = CM.series!(primVals,transpose(dropdims(mean(meanResults[zParam],dims=3),dims=3)),
            labels = seriesLabels, linewidth = 1, color = plotColors);
        CM.vspan!(1,3,color=unboundedColor,
            label="species richness\nunbounded",labelfont = :bold,overdraw = true)
        
        CM.rowgap!(gridLayout,10)
        CM.linkxaxes!(axb,axt)
    
        for (label, layout) in zip(["A", "B"], [fig[1,1], fig[2,1]])
            CM.Label(layout[1, 1, CM.TopLeft()], label,
                #fontsize = 26,
                font = :bold,
                padding = (0, 25, 5, 0),
                halign = :left)
        end

        save("data/finalPlots/primExtEatsNeeds.$(plotDataType)", fig, pt_per_unit = 1);
    else
        error("Plotting library \"$plotLib\" not supported.")
    end

    #meanPlot = plot(primVals[10:end], dropdims(mean(meanResults[:specRich],dims=3),dims=3)[10:end,:],
    #ylabel = zParamLongNames[:specRich],
    ##title = "Average mean $(zParamLongNames[zParam]) in itterations $(maxits - 1500) to $maxits",
    #yaxis = :log, labels = "specRich: " .* seriesLabels, legend = :topleft)
#
    #meanPlot = plot!(twinx(),primVals[10:end], dropdims(mean(meanResults[:meanEats],dims=3),dims=3)[10:end,:],
    #xlabel = "primary extinction rate", ylabel = zParamLongNames[:meanEats],
    ##title = "Average mean $(zParamLongNames[zParam]) in itterations $(maxits - 1500) to $maxits",
    #yaxis = :log, labels = "meanEats: " .* seriesLabels,c = [:green :black])
#
    #Plots.savefig(meanPlot,"data/$simulationName/plots/mean_specRichAndMeanEatsPlot.html");
end