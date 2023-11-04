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
@everywhere using GLM
@everywhere using DataFrames

@everywhere @enum NetworkGrowthType begin
    uncertain = 0
    boundedUncertain = 1
    boundedStable = 2
    boundedUnstable = 3
    unboundedUncertain = 4
    undboundedLinear = 5
    unboundedExponential = 6
    trivial = 7
end

isBounded(netGrowthType::NetworkGrowthType) =
    netGrowthType in (trivial,boundedUncertain,boundedStable, boundedUnstable)

isUnbounded(netGrowthType::NetworkGrowthType) =
    netGrowthType in (unboundedUncertain, undboundedLinear, unboundedExponential)



function boundednessHeatmap()
    include("src/set_up_params.jl")

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

    bounded = SharedArray{Bool}((numParams,numParams,nRepets))
    netGrowthTypes = SharedArray{NetworkGrowthType}((numParams,numParams,nRepets))

    maxSlope = 0.1
    maxRecursionDepth = 5

    @everywhere function checkBoundedness(sd,rates,bounded,netGrowthTypes,primInd,rPrimExt,secInd,rSecExt,repet,maxSlope,maxRecursionDepth,
        maxits,cm,cn,ce,cPred,diverse,restrict_colonization,recursionDepth=0;offset = 1)
        linReg = lm(@formula(specRich ~ clock), DataFrame(clock = sd.clock[offset:end], specRich = sd.specRich[offset:end]))

        coefTable = GLM.coeftable(linReg).cols

        lower95 = coefTable[5][2]
        upper95 = coefTable[6][2]

        if lower95 > maxSlope
            bounded[primInd,secInd,repet] = false
            netGrowthTypes[primInd,secInd,repet] = unboundedUncertain
        elseif upper95 < maxSlope# && lower95 > -maxSlope
            bounded[primInd,secInd,repet] = true
            netGrowthTypes[primInd,secInd,repet] = boundedUncertain
        else
            if recursionDepth < maxRecursionDepth
                sd,_ = assemblyevo(sd.poolnet, rates, maxits, cn,cm,ce,cpred, diverse, restrict_colonization, sd.colnet)
                checkBoundedness(sd,rates,bounded,netGrowthTypes,primInd,rPrimExt,secInd,rSecExt,repet,maxSlope,maxRecursionDepth,
                    maxits,cm,cn,ce,cPred,diverse,restrict_colonization,recursionDepth+1)
            else
                netGrowthTypes[primInd,secInd,repet] = uncertain
            end
        end
    end

    loopVars = [(primInd,rPrimExt,secInd,rSecExt,repetition) for (primInd,rPrimExt) in enumerate(paramVals)
        for (secInd,rSecExt) in enumerate(paramVals) for repetition in repets];

    @sync @distributed for (primInd,rPrimExt,secInd,rSecExt,repetition) in loopVars
        sd,rates0 = load("data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rSecExt))_repet=$repetition.jld2", "simulationData","rates0")

        checkBoundedness(sd,rates0,bounded,netGrowthTypes,primInd,rPrimExt,secInd,rSecExt,repetition,maxSlope,maxRecursionDepth,
            maxits,cm,cn,ce,cpred,diverse,restrict_colonization;offset = 2_000)
        #nets[primInd,secInd] = sd
        #plots[primInd,secInd] = plot(sd.clock,sd.specRich,title="(rPrimExt,rSecExt) = $((rPrimExt,rSecExt))",ylabel = "specRich", xlabel = "clock", show = false)
    end

    plotlyjs()

    boundedPlot = Plots.heatmap(string.(paramVals), string.(paramVals),
        dropdims(mapslices(page->count(isBounded,page), netGrowthTypes, dims = 3),dims=3),
        size = (720,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
        title = "Number of runs categorized as bounded.", xticks = :all,
        yticks = :all, xrotation = 60, show = false)

    Plots.savefig(boundedPlot,"data/$simulationName/plots/boundedPlot.html");

    unboundedPlot = Plots.heatmap(string.(paramVals), string.(paramVals),
        dropdims(mapslices(page->count(isUnbounded,page), netGrowthTypes, dims = 3),dims=3),
        size = (720,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
        title = "Number of runs categorized as unbounded.", xticks = :all,
        yticks = :all, xrotation = 60, show = false)

    Plots.savefig(unboundedPlot,"data/$simulationName/plots/unboundedPlot.html");

    uncertainPlot = Plots.heatmap(string.(paramVals), string.(paramVals),
        dropdims(mapslices(page->count(==(uncertain),page), netGrowthTypes, dims = 3),dims=3),
        size = (720,720), xlabel = "secondary extinction rate", ylabel = "primary extinction rate",
        title = "Number of runs that could not be categorized.", xticks = :all,
        yticks = :all, xrotation = 60, show = false)

    Plots.savefig(uncertainPlot,"data/$simulationName/plots/uncertainPlot.html");
end

function boundednessPrimExtLinePlot(onlyPlots=false)
    include("src/set_up_params.jl")
    maxits = 30_000
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

    if !onlyPlots
        bounded = SharedArray{Bool}((nPrimVals,nSecVals,nRepets))
        netGrowthTypes = SharedArray{NetworkGrowthType}((nPrimVals,nSecVals,nRepets))

        maxSlope = 0.1
        maxRecursionDepth = 0

        @everywhere function checkBoundedness(sd,rates,bounded,netGrowthTypes,rPrimInd,secInd,repet,maxSlope,maxRecursionDepth,
            maxits,cm,cn,ce,cPred,diverse,restrict_colonization,recursionDepth=0;offset = 1)
            linReg = lm(@formula(specRich ~ clock), DataFrame(clock = sd.clock[offset:end], specRich = sd.specRich[offset:end]))

            coefTable = GLM.coeftable(linReg).cols

            lower95 = coefTable[5][2]
            upper95 = coefTable[6][2]

            if lower95 > maxSlope
                bounded[rPrimInd,secInd,repet] = false
                netGrowthTypes[rPrimInd,secInd,repet] = unboundedUncertain
            elseif upper95 < maxSlope# && lower95 > -maxSlope
                bounded[rPrimInd,secInd,repet] = true
                netGrowthTypes[rPrimInd,secInd,repet] = boundedUncertain
            else
                if recursionDepth < maxRecursionDepth
                    sd,_ = assemblyevo(sd.poolnet, rates, maxits, cn,cm,ce,cpred, diverse, restrict_colonization, sd.colnet)
                    checkBoundedness(sd,rates,bounded,netGrowthTypes,rPrimInd,secInd,repet,maxSlope,maxRecursionDepth,
                        maxits,cm,cn,ce,cPred,diverse,restrict_colonization,recursionDepth+1)
                else
                    netGrowthTypes[rPrimInd,secInd,repet] = uncertain
                end
            end
        end

        loop_vars = [(primInd,rPrimExt,secInd,rSecExt,repetition) for (primInd,rPrimExt) in enumerate(primVals)
        for (secInd,rSecExt) in enumerate(secVals) for repetition in repets];

        @sync @distributed for (rPrimInd,rPrimExt,secInd,rSecExt,repetition) in loop_vars
            fileName = "data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rSecExt))_repet=$repetition.jld2" 
            sd,rates0 = load(fileName,"simulationData","rates0")

            checkBoundedness(sd,rates0,bounded,netGrowthTypes,rPrimInd,secInd,repetition,maxSlope,maxRecursionDepth,
                maxits,cm,cn,ce,cpred,diverse,restrict_colonization;offset = 15_000)
        end
        jldsave("data/$(simulationName)/BoundedResults.jld2",compress; bounded, netGrowthTypes)
    else
        bounded, netGrowthTypes = load("data/$(simulationName)/BoundedResults.jld2","bounded","netGrowthTypes")
    end


    plotlyjs()

    seriesLabels = reshape(["rSecExt = $rSecExt" for rSecExt in secVals], (1,nSecVals))

    boundedPlot = Plots.plot(primVals, dropdims(mapslices(page->count(isBounded,page), netGrowthTypes, dims = 3),dims=3),
        xlabel = "rPrim", ylabel = "number of runs categorized as bounded",
        title = "Number of runs categorized as bounded.",
        labels = seriesLabels, show = false)

    Plots.savefig(boundedPlot,"data/$simulationName/plots/boundedPlot.html");

    unboundedPlot = Plots.plot(primVals, dropdims(mapslices(page->count(isUnbounded,page), netGrowthTypes, dims = 3),dims=3),
        xlabel = "rPrim", ylabel = "number of runs categorized as unbounded",
        title = "Number of runs categorized as unbounded.",
        labels = seriesLabels, show = false)

    Plots.savefig(unboundedPlot,"data/$simulationName/plots/unboundedPlot.html");

    uncertainPlot = Plots.plot(primVals, dropdims(mapslices(page->count(==(uncertain),page), netGrowthTypes, dims = 3),dims=3),
        xlabel = "rPrim", ylabel = "number of runs categorized as uncertain",
        title = "Number of runs categorized as uncertain.",
        labels = seriesLabels, show = false)

    Plots.savefig(uncertainPlot,"data/$simulationName/plots/uncertainPlot.html");
end

function boundednessEvoLinePlot()
    include("src/set_up_params.jl")
    maxits = 30_000
    #prepare everything for a simulation consisting of the variation of a parmeter
    paramName = "rPrimExt,rSecExt,evoInd"
    simulationName = "rEvoLinePlotconstGlobExt"#"rEvoInDepth";        #specify the name of the simulation
    mkpath("data/$(simulationName)/runs");    #make a folder with that name in the Data folder including plots subfolder
    mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder including plots subfolder
    compress::Bool = true;  #should the data be compressed before storing?
 
    primVals = [3.5]       #specify the parameter values that shall be simulated
    nPrimVals = length(primVals)
    secVals = [1.,10.]
    nSecVals = length(secVals)
    evoVals = exp10.(range(log10(0.005), stop=log10(2), length=40))[10:30]
    nEvoVals = length(evoVals)
    evoInds = [6,28,37];    #hand picked indices of values of special interest 
    nEvoInds = length(evoInds)
    nRepets = 50
    repets = 1:nRepets

    bounded = SharedArray{Bool}((nEvoVals,nSecVals,nRepets))
    netGrowthTypes = SharedArray{NetworkGrowthType}((nEvoVals,nSecVals,nRepets))

    maxSlope = 0.1
    maxRecursionDepth = 0

    @everywhere function checkBoundedness(sd,rates,bounded,netGrowthTypes,rEvoInd,secInd,repet,maxSlope,maxRecursionDepth,
        maxits,cm,cn,ce,cPred,diverse,restrict_colonization,recursionDepth=0;offset = 1)
        linReg = lm(@formula(specRich ~ clock), DataFrame(clock = sd.clock[offset:end], specRich = sd.specRich[offset:end]))

        coefTable = GLM.coeftable(linReg).cols

        lower95 = coefTable[5][2]
        upper95 = coefTable[6][2]

        if lower95 > maxSlope
            bounded[rEvoInd,secInd,repet] = false
            netGrowthTypes[rEvoInd,secInd,repet] = unboundedUncertain
        elseif upper95 < maxSlope# && lower95 > -maxSlope
            bounded[rEvoInd,secInd,repet] = true
            netGrowthTypes[rEvoInd,secInd,repet] = boundedUncertain
        else
            if recursionDepth < maxRecursionDepth
                sd,_ = assemblyevo(sd.poolnet, rates, maxits, cn,cm,ce,cpred, diverse, restrict_colonization, sd.colnet)
                checkBoundedness(sd,rates,bounded,netGrowthTypes,rEvoInd,secInd,repet,maxSlope,maxRecursionDepth,
                    maxits,cm,cn,ce,cPred,diverse,restrict_colonization,recursionDepth+1)
            else
                netGrowthTypes[rEvoInd,secInd,repet] = uncertain
            end
        end
    end


    loop_vars = [(rEvoInd,rPrimExt,secInd,rSecExt,repetition) for repetition in repets
        for rPrimExt in primVals for (secInd,rSecExt) in enumerate(secVals) for rEvoInd in eachindex(evoVals)];

    @sync @distributed for (rEvoInd,rPrimExt,secInd,rSecExt,repetition) in loop_vars
        sd,rates0 = load("data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rSecExt),$(rEvoInd))_repet=$repetition.jld2",
            "simulationData","rates0")

        checkBoundedness(sd,rates0,bounded,netGrowthTypes,rEvoInd,secInd,repetition,maxSlope,maxRecursionDepth,
            maxits,cm,cn,ce,cpred,diverse,restrict_colonization;offset = 15_000)
        #nets[primInd,rEvoInd] = sd
        #plots[primInd,rEvoInd] = plot(sd.clock,sd.specRich,title="(rPrimExt,rEvo) = $((rPrimExt,rEvo))",ylabel = "specRich", xlabel = "clock", show = false)
    end

    plotlyjs()

    seriesLabels = reshape(["rSecExt = $rSecExt" for rSecExt in secVals], (1,nSecVals))

    boundedPlot = Plots.plot(evoVals, dropdims(mapslices(page->count(isBounded,page), netGrowthTypes, dims = 3),dims=3),
        xlabel = "rEvo", ylabel = "number of runs categorized as bounded",
        title = "Number of runs categorized as bounded.", xaxis = :log,
        labels = seriesLabels, show = false)

    Plots.savefig(boundedPlot,"data/$simulationName/plots/boundedPlot.html");

    unboundedPlot = Plots.plot(evoVals, dropdims(mapslices(page->count(isUnbounded,page), netGrowthTypes, dims = 3),dims=3),
        xlabel = "rEvo", ylabel = "number of runs categorized as unbounded",
        title = "Number of runs categorized as unbounded.", xaxis = :log,
        labels = seriesLabels, show = false)

    Plots.savefig(unboundedPlot,"data/$simulationName/plots/unboundedPlot.html");

    uncertainPlot = Plots.plot(evoVals, dropdims(mapslices(page->count(==(uncertain),page), netGrowthTypes, dims = 3),dims=3),
        xlabel = "rEvo", ylabel = "number of runs categorized as uncertain",
        title = "Number of runs categorized as uncertain.", xaxis = :log,
        labels = seriesLabels, show = false)

    Plots.savefig(uncertainPlot,"data/$simulationName/plots/uncertainPlot.html");
end

function boundednessHeatmapEvo(;couplingFactor = .5)
    include("src/set_up_params.jl")

    #prepare everything for a simulation consisting of the variation of a parmeter
    paramName = "(rPrimExt,rEvo)"
    simulationName = "heatmapPrimExtVs(REvoCoupledRGlob)";        #specify the name of the simulation
    mkpath("data/$(simulationName)/plots");     #make a folder with that name in the Data folder including plots subfolder
    mkpath("data/$(simulationName)/runs"); 
    compress::Bool = true;  #should the data be compressed before storing?
    
    primVals = [1,1.5,2,2.5,3,4,5,10,15,20]     #specify the parameter values that shall be simulated
    rEvoVals = [0.005, 0.01, 0.016, 0.022, 0.03, 0.035, 0.045, 0.06, 0.08, 0.22, 0.3, 0.35, 0.45, 0.55, 0.7, 1.0, 2.0]  #0.015:.0005:0.022
    numPrimVals = length(primVals)
    numREvoVals = length(rEvoVals)
    nRepets = 50
    repets = 1:nRepets

    bounded = SharedArray{Bool}((numPrimVals,numREvoVals,nRepets))
    netGrowthTypes = SharedArray{NetworkGrowthType}((numPrimVals,numREvoVals,nRepets))

    maxSlope = 0.1
    maxRecursionDepth = 0

    @everywhere function checkBoundedness(sd,rates,bounded,netGrowthTypes,primInd,rPrimExt,rEvoInd,rEvo,repet,maxSlope,maxRecursionDepth,
        maxits,cm,cn,ce,cPred,diverse,restrict_colonization,recursionDepth=0;offset = 1)
        linReg = lm(@formula(specRich ~ clock), DataFrame(clock = sd.clock[offset:end], specRich = sd.specRich[offset:end]))

        coefTable = GLM.coeftable(linReg).cols

        lower95 = coefTable[5][2]
        upper95 = coefTable[6][2]

        if lower95 > maxSlope
            bounded[primInd,rEvoInd,repet] = false
            netGrowthTypes[primInd,rEvoInd,repet] = unboundedUncertain
        elseif upper95 < maxSlope# && lower95 > -maxSlope
            bounded[primInd,rEvoInd,repet] = true
            netGrowthTypes[primInd,rEvoInd,repet] = boundedUncertain
        else
            if recursionDepth < maxRecursionDepth
                sd,_ = assemblyevo(sd.poolnet, rates, maxits, cn,cm,ce,cpred, diverse, restrict_colonization, sd.colnet)
                checkBoundedness(sd,rates,bounded,netGrowthTypes,primInd,rPrimExt,rEvoInd,rEvo,repet,maxSlope,maxRecursionDepth,
                    maxits,cm,cn,ce,cPred,diverse,restrict_colonization,recursionDepth+1)
            else
                netGrowthTypes[primInd,rEvoInd,repet] = uncertain
            end
        end
    end


    loop_vars = [(primInd,rPrimExt,rEvoInd,rEvo,repetition) for repetition in repets
        for (primInd,rPrimExt) in enumerate(primVals) for (rEvoInd,rEvo) in enumerate(rEvoVals)];

    @sync @distributed for (primInd,rPrimExt,rEvoInd,rEvo,repetition) in loop_vars
        sd,rates0 = load("data/$(simulationName)/runs/$(paramName)=($(rPrimExt),$(rEvo))_repet=$repetition.jld2",
            "simulationData","rates0")

        checkBoundedness(sd,rates0,bounded,netGrowthTypes,primInd,rPrimExt,rEvoInd,rEvo,repetition,maxSlope,maxRecursionDepth,
            maxits,cm,cn,ce,cpred,diverse,restrict_colonization;offset = 2_000)
        #nets[primInd,rEvoInd] = sd
        #plots[primInd,rEvoInd] = plot(sd.clock,sd.specRich,title="(rPrimExt,rEvo) = $((rPrimExt,rEvo))",ylabel = "specRich", xlabel = "clock", show = false)
    end

    plotlyjs()

    boundedPlot = Plots.heatmap("(".*string.(rEvoVals).*", ".*string.(rEvoVals.*couplingFactor).*")",
        string.(primVals), dropdims(mapslices(page->count(isBounded,page), netGrowthTypes, dims = 3),dims=3),
        size = (720,720), xlabel = "(rEvo,rExt)", ylabel = "primary extinction rate",
        title = "Number of runs categorized as bounded.", xticks = :all,
        yticks = :all, xrotation = 60, show = false)

    Plots.savefig(boundedPlot,"data/$simulationName/plots/boundedPlot.html");

    unboundedPlot = Plots.heatmap("(".*string.(rEvoVals).*", ".*string.(rEvoVals.*couplingFactor).*")",
    string.(primVals), dropdims(mapslices(page->count(isUnbounded,page), netGrowthTypes, dims = 3),dims=3),
    size = (720,720), xlabel = "(rEvo,rExt)", ylabel = "primary extinction rate",
    title = "Number of runs categorized as unbounded.", xticks = :all,
    yticks = :all, xrotation = 60, show = false)

    Plots.savefig(unboundedPlot,"data/$simulationName/plots/unboundedPlot.html");

    uncertainPlot = Plots.heatmap("(".*string.(rEvoVals).*", ".*string.(rEvoVals.*couplingFactor).*")",
    string.(primVals), dropdims(mapslices(page->count(==(uncertain),page), netGrowthTypes, dims = 3),dims=3),
    size = (720,720), xlabel = "(rEvo,rExt)", ylabel = "primary extinction rate",
    title = "Number of runs that could not be categorized.", xticks = :all,
    yticks = :all, xrotation = 60, show = false)

    Plots.savefig(uncertainPlot,"data/$simulationName/plots/uncertainPlot.html");
end