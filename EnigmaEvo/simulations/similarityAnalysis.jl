if !@isdefined(NO_STANDALONE) || NO_STANDALONE == false #for usage in scripts where only the function is needed and the rest will already be loaded
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
end

""""
    doSimilarityAnalysis(onlyPlots=false,fromResults=false;couplingFactor=.5, plotLib=:Plots)    

    Performs an analysis of the pairwise similarities of the species in a set of simulations measured with a modified Sørensen-Dice coefficient. Creates line plots of results.

    # Arguments
    'onlyPlots::Bool': If true assumes the simulations to be finished and only performs data analysis from files and plots. 
    'fromResults::Bool': If true loads final results and plots. Overrides 'onlyPlots'.
    'couplingFactor::Float': Defines the coupling factor of evolutionary and global extinction rate.'
    'plotLib::Symbol': Symbolic name of plotting library. Supported are :Plots and :Makie
"""
function doSimilarityAnalysis(onlyPlots=false,fromResults=false;couplingFactor=.5, plotLib=:Plots)   #supposedly its better to wrap stuff in functions and not use global variables if possible(not sure if typed ones are a problem)
    if fromResults && !onlyPlots
        error("If only results file should be used it implies that only plots can be created.")
    end

    include("src/set_up_params.jl");
    maxits = 30_000
    #prepare everything for a simulation consisting of the variation of a parmeter
    paramName = "rPrimExt,rSecExt,evoInd"
    simulationName = "rEvoLinePlot";        #specify the name of the simulation
    mkpath("data/$(simulationName)/runs");    #make a folder with that name in the Data folder including plots subfolder
    mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder including plots subfolder
    compress::Bool = true;  #should the data be compressed before storing?


    primVals = [3.5]       #specify the parameter values that shall be simulated
    nPrimVals = length(primVals)
    secVals = [1.,10.]
    nSecVals = length(secVals)
    evoVals = exp10.(range(log10(0.005), stop=log10(2), length=40))[1:37]
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
        :SDsPool => "Species in pool",
        :SDsCol => "Species in colony",
        :SDsColInPool => "Species in colony as subset of pool",
        :SDsPoolWithoutCol => "Species in pool without species in colony",
        :SDsPoolOnlyEats => "Species in pool reduced to eats",
        :SDsColOnlyEats => "Species in colony reduced to eats",
        :SDsColInPoolOnlyEats => "Species in colony as subset of pool reduced to eats",
        :SDsPoolWithoutColOnlyEats => "Species in pool without species in colony reduced to eats"
    )

    if !fromResults
        SDmeans = Dict{Symbol, SharedArray{Float64}}()                                                                                                   
        SDstds = Dict{Symbol, SharedArray{Float64}}()                                                                                                                                          
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
                continue
            end
            if sd.runFinished    #did the run end prematurely?
                SDsPool,specIdsPool = modifiedSDSimilarities(sd.poolnet);
                SDsCol,specIdsCol = modifiedSDSimilarities(sd.colnet);
                SDsColInPool,specIdsColInPool = modifiedSDSimilarities(sd.poolnet,collect(sd.colnet.spec));
                SDsPoolWithoutCol,specIdsPoolWithoutCol = modifiedSDSimilarities(sd.poolnet,collect(setdiff(sd.poolnet.spec,sd.colnet.spec)));
                SDsPoolOnlyEats,specIdsPoolOnlyEats = modifiedSDSimilarities(sd.poolnet,:eat);
                SDsColOnlyEats,specIdsColOnlyEats = modifiedSDSimilarities(sd.colnet,:eat);
                SDsColInPoolOnlyEats,specIdsColInPoolOnlyEats = modifiedSDSimilarities(sd.poolnet,collect(sd.colnet.spec),:eat);
                SDsPoolWithoutColOnlyEats,specIdsPoolWithoutColOnlyEats = modifiedSDSimilarities(sd.poolnet,collect(setdiff(sd.poolnet.spec,sd.colnet.spec)),:eat);
                
                SDmeans[:SDsPool][primInd,evoInd,secInd,repetition] = nanmean(SDsPool)
                SDmeans[:SDsPoolOnlyEats][primInd,evoInd,secInd,repetition] = nanmean(SDsPoolOnlyEats)
                SDmeans[:SDsCol][primInd,evoInd,secInd,repetition] = nanmean(SDsCol)
                SDmeans[:SDsColOnlyEats][primInd,evoInd,secInd,repetition] = nanmean(SDsColOnlyEats)
                SDmeans[:SDsColInPool][primInd,evoInd,secInd,repetition] = nanmean(SDsColInPool)
                SDmeans[:SDsColInPoolOnlyEats][primInd,evoInd,secInd,repetition] = nanmean(SDsColInPoolOnlyEats)
                SDmeans[:SDsPoolWithoutCol][primInd,evoInd,secInd,repetition] = nanmean(SDsPoolWithoutCol)
                SDmeans[:SDsPoolWithoutColOnlyEats][primInd,evoInd,secInd,repetition] = nanmean(SDsPoolWithoutColOnlyEats)
                SDstds[:SDsPool][primInd,evoInd,secInd,repetition] = nanstd(SDsPool)
                SDstds[:SDsPoolOnlyEats][primInd,evoInd,secInd,repetition] = nanstd(SDsPoolOnlyEats)
                SDstds[:SDsCol][primInd,evoInd,secInd,repetition] = nanstd(SDsCol)
                SDstds[:SDsColOnlyEats][primInd,evoInd,secInd,repetition] = nanstd(SDsColOnlyEats)
                SDstds[:SDsColInPool][primInd,evoInd,secInd,repetition] = nanstd(SDsColInPool)
                SDstds[:SDsColInPoolOnlyEats][primInd,evoInd,secInd,repetition] = nanstd(SDsColInPoolOnlyEats)
                SDstds[:SDsPoolWithoutCol][primInd,evoInd,secInd,repetition] = nanstd(SDsPoolWithoutCol)
                SDstds[:SDsPoolWithoutColOnlyEats][primInd,evoInd,secInd,repetition] = nanstd(SDsPoolWithoutColOnlyEats)
            end
        end
        jldsave("data/$(simulationName)/similarityResults.jld2", compress; SDmeans, SDstds)
    else
        SDmeans, SDstds = 
            load("data/$(simulationName)/similarityResults.jld2", "SDmeans", "SDstds")
    end

    seriesLabels = reshape([L"r_{2\degree} = %$rSecExt" for rSecExt in secVals], (1,nSecVals))
    xLabel = L"r_{evo}=2r_g"
    yLabelMean = "mean Sørensen-Dice coefficient"
    yLabelStd = "standard deviation of Sørensen-Dice coefficient}"
    xScale = log10
    if plotLib == :Plots
        plotlyjs()
        xLabel = "evolutionary rate ($(1/couplingFactor)*global extinction rate)"
        seriesLabels = reshape(["rSecExt = $rSecExt" for rSecExt in secVals], (1,nSecVals))
        for (primInd,rPrimExt) in enumerate(primVals)
            plotPath = "data/$simulationName/plots/rPrimExt=$(rPrimExt)"
            mkpath(plotPath)

            for category in categories
                meanPlot = plot(evoVals, dropdims(mean(selectdim(SDmeans[category],1,primInd),dims=3),dims=3),
                    xlabel = xLabel, ylabel = yLabelMean,
                    title = "Mean Sorensen-Dice Similarity Coefficient for all $(categoryLongNames[category]).",
                    xaxis = :log, labels = seriesLabels, legend = :topright)

                Plots.savefig(meanPlot,"$(plotPath)/meanSD_$(String(category))Plot.svg");

                if any(isnan,SDmeans[category])
                    meanPlot = plot(evoVals, NaNStatistics.nanmean(selectdim(SDmeans[category],1,primInd),dim=3),
                    xlabel = xLabel, ylabel = yLabelMean,
                    title = "Mean Sorensen-Dice Similarity Coefficient for all $(categoryLongNames[category]) - ignoring NaNs.", 
                    xaxis = :log, labels = seriesLabels, legend = :topleft)

                    Plots.savefig(meanPlot,"$(plotPath)/meanSD_$(String(category))Plot_ignoreNaNs.svg");
                end

                stdPlot = plot(evoVals, dropdims(mean(selectdim(SDmeans[category],1,primInd),dims=3),dims=3),
                    xlabel = xLabel, ylabel = yLabelStd,
                    title = "Standard deviation of Sorensen-Dice Similarity Coefficient for all $(categoryLongNames[category]).",
                    xaxis = :log, labels = seriesLabels, legend = :topright)

                Plots.savefig(stdPlot,"$(plotPath)/stdSD_$(String(category))Plot.svg");

                if any(isnan,SDstds[category])
                    stdPlot = plot(evoVals, NaNStatistics.nanmean(selectdim(SDmeans[category],1,primInd),dim=3),
                        xlabel = xLabel, ylabel = yLabelStd,
                        title = "Standard deviation of Sorensen-Dice Similarity Coefficient for all $(categoryLongNames[category]) - ignoring NaNs.",  
                        xaxis = :log, labels = seriesLabels, legend = :topleft)

                    Plots.savefig(stdPlot,"$(plotPath)/stdSD_$(String(category))Plot_ignoreNaNs.svg");
                end
            end
        end
    elseif plotLib == :Makie
        unboundedMin = evoVals[33]
        unboundedMax = evoVals[37]
        axisLimits = ((nothing,unboundedMax),nothing)
        for (primInd,rPrimExt) in enumerate(primVals)
            plotPath = "data/$simulationName/plots/rPrimExt=$(rPrimExt)"
            mkpath(plotPath)

            for category in categories
                fig = CM.Figure(resolution = plotSizePt, fontsize = fontSize, figure_padding = figurePadding);
                ax = CM.Axis(fig[1,1];xlabel = xLabel, ylabel = yLabelMean,
                            xscale = xScale, yscale = identity, limits = axisLimits);
                se = CM.series!(ax,evoVals,
                    transpose(dropdims(mean(selectdim(SDmeans[category],1,primInd),dims=3),dims=3)),
                    labels = seriesLabels, linewidth = lineWidth, color = plotColors);
                CM.vspan!(unboundedMin,unboundedMax,color=unboundedColor,
                    label="species richness\nunbounded",labelfont = :bold,overdraw = true)
                CM.axislegend(ax,patchsize=legendPatchSize, position = legendPosition)   
        
                save("$(plotPath)/meanSD_$(String(category))Plot.$(plotDataType)", fig, pt_per_unit = 1);

                if any(isnan,SDmeans[category])
                    fig = CM.Figure(resolution = plotSizePt, fontsize = fontSize, figure_padding = figurePadding);
                    ax = CM.Axis(fig[1,1];xlabel = xLabel, ylabel = yLabelMean,
                                xscale = xScale, yscale = identity, limits = axisLimits);
                    se = CM.series!(ax,evoVals,
                        transpose(NaNStatistics.nanmean(selectdim(SDmeans[category],1,primInd),dim=3)),
                        labels = seriesLabels, linewidth = lineWidth, color = plotColors);
                    CM.vspan!(unboundedMin,unboundedMax,color=unboundedColor,
                        label="species richness\nunbounded",labelfont = :bold,overdraw = true)
                    CM.axislegend(ax,patchsize=legendPatchSize, position = legendPosition)   
            
                    save("$(plotPath)/meanSD_$(String(category))Plot_ignoreNaNs.$(plotDataType)", fig, pt_per_unit = 1);
                end

                fig = CM.Figure(resolution = plotSizePt, fontsize = fontSize, figure_padding = figurePadding);
                ax = CM.Axis(fig[1,1];xlabel = xLabel, ylabel = yLabelStd,
                            xscale = xScale, yscale = identity, limits = axisLimits);
                se = CM.series!(ax,evoVals,
                    transpose(dropdims(mean(selectdim(SDmeans[category],1,primInd),dims=3),dims=3)),
                    labels = seriesLabels, linewidth = lineWidth, color = plotColors);
                CM.vspan!(unboundedMin,unboundedMax,color=unboundedColor,
                    label="species richness\nunbounded",labelfont = :bold,overdraw = true)
                CM.axislegend(ax,patchsize=legendPatchSize, position = legendPosition)   
        
                save("$(plotPath)/stdSD_$(String(category))Plot.$(plotDataType)", fig, pt_per_unit = 1);

                if any(isnan,SDstds[category])
                    fig = CM.Figure(resolution = plotSizePt, fontsize = fontSize, figure_padding = figurePadding);
                    ax = CM.Axis(fig[1,1];xlabel = xLabel, ylabel = yLabelStd,
                                xscale = xScale, yscale = identity, limits = axisLimits);
                    se = CM.series!(ax,evoVals,
                        transpose(NaNStatistics.nanmean(selectdim(SDmeans[category],1,primInd),dim=3)),
                        labels = seriesLabels, linewidth = lineWidth, color = plotColors);
                    CM.vspan!(unboundedMin,unboundedMax,color=unboundedColor,
                        label="species richness\nunbounded",labelfont = :bold,overdraw = true)
                    CM.axislegend(ax,patchsize=legendPatchSize, position = legendPosition)   
            
                    save("$(plotPath)/stdSD_$(String(category))Plot_ignoreNaNs.$(plotDataType)", fig, pt_per_unit = 1);
                end
            end

            #final plots
            fig = CM.Figure(resolution = triplePlotSizePt, fontsize = fontSize, figure_padding = figurePadding);
            gridLayout = fig.layout
            category = :SDsCol
            axt = CM.Axis(gridLayout[1,1];ylabel = yLabelMean, ylabelpadding = 17,
                            xscale= log10, yscale = identity, limits = axisLimits);
            CM.hidexdecorations!(axt,grid = false)
            se = CM.series!(axt,evoVals,
                transpose(NaNStatistics.nanmean(selectdim(SDmeans[category],1,primInd),dim=3)),
                labels = seriesLabels, linewidth = lineWidth, color = plotColors);
            CM.vspan!(unboundedMin, unboundedMax, color=unboundedColor,
                label="species richness\nunbounded", labelfont = :bold, overdraw = true)
            CM.axislegend(axt,patchsize=legendPatchSize, position = legendPosition)   

            axb = CM.Axis(gridLayout[2,1];xlabel = xLabel,ylabel = "Pearson correlation coefficient of\nDSC and species richness",
                xscale = log10, yscale = identity, limits = axisLimits);
            se = CM.series!(axb,evoVals, transpose([corSec1 corSec10]),
                labels = seriesLabels, linewidth = lineWidth, color = plotColors)
            CM.vspan!(unboundedMin,unboundedMax,color=unboundedColor,
                label="species richness\nunbounded",labelfont = :bold,overdraw = true)
            
            CM.linkxaxes!(axb,axt)

            evoInd = 12
            ax3 = CM.Axis(gridLayout[3,1];ylabel = "mean Sørensen-Dice coefficient",
                xlabel = "species richness",ylabelpadding = 12.5);
            CM.scatter!(ax3,meanResults[:specRich][1,13,1,:],SDmeans[:SDsCol][1,12,1,:],
                colormap = plotColors,markersize = 5, label = L"r_{2\degree} = 1.0")
            CM.scatter!(ax3,meanResults[:specRich][1,13,2,:],SDmeans[:SDsCol][1,12,2,:],
                colormap = plotColors,markersize = 5, label = L"r_{2\degree} = 10.0")
            CM.axislegend(patchsize=legendPatchSize)

            CM.rowgap!(gridLayout,10)

            for (label, layout) in zip(["A", "B", "C"], [fig[1,1], fig[2,1], fig[3,1]])
                CM.Label(layout[1, 1, CM.TopLeft()], label,
                    #fontsize = 26,
                    font = :bold,
                    padding = (0, 32, 5, 0),
                    halign = :left)
            end
            save("data/finalPlots/evoSimilarity.$(plotDataType)", fig, pt_per_unit = 1);
        end
    else
        error("Plotting library \"$plotLib\" not supported.")
    end
end