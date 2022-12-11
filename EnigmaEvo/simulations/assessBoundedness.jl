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
using GLM
using DataFrames

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

    const maxSlope = 0.1
    const maxRecursionDepth = 5

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

    loopVars = [(primInd,rPrimExt,secInd,rSecExt,repetition) for (primInd,rPrimExt) in enumerate(paramVals)
        for (secInd,rSecExt) in enumerate(paramVals) for repetition in repets];

    @sync @distributed for (primInd,rPrimExt,secInd,rSecExt,repetition) in loopVars
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
end