if homedir() == "/home/z840"    #for downward compatibility ;)
    localPath::String = "$(homedir())/2019_Lego_Evo/EnigmaEvo/src/";
elseif isfile("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl")
    localPath::String = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/";
else                            #othserwise use relative paths
    localPath::String = "src/";
end;

using Distributed
addprocs(80)

include(localPath*"loadfuncs.jl");

using SharedArrays

function simulation()   #supposedly its better to wrap stuff in functions and not use global variables if possible(not sure if typed ones are a problem)
    include("src/set_up_params.jl");
    maxits = 10_000
    #prepare everything for a simulation consisting of the variation of a parmeter
    paramName = "(rPrimExt,rSecExt)"
    simulationName = "varyExtinctionsForTrophLevel";        #specify the name of the simulation
    mkpath("data/$(simulationName)/plots");    #make a folder with that name in the Data folder including plots subfolder
    compress::Bool = true;  #should the data be compressed before storing?

    jldsave("data/$(simulationName)/parameters.jld2",compress;
        rates0, maxits, cm,ce,cpred, diverse, restrict_colonization,
        logging,S,lambda,SSprobs,SOprobs,nBasalRes)

    paramVals = (1:.5:10).^2;       #specify the parameter values that shall be simulated
    numParams = length(paramVals)
    nRepets = 50
    repets = 1:nRepets
               
    loop_vars = [(primInd,rPrimExt,secInd,rSecExt,repetition) for (primInd,rPrimExt) in enumerate(paramVals)
        for (secInd,rSecExt) in enumerate(paramVals) for repetition in repets];

    heatMapMax = SharedArray{Float64}((numParams,numParams,nRepets));
    heatMapMean = SharedArray{Float64}((numParams,numParams,nRepets));
    @sync @distributed for (primInd,rPrimExt,secInd,rSecExt,repetition) in loop_vars
        initpoolnet::ENIgMaGraph = setUpPool(S,lambda,nBasalRes,SSprobs,SOprobs,diverse);
        rates0 = (rc = 1., rprimext = rPrimExt, rsecext = rSecExt, reo = 1., revo = 0.035, rext = 0.016);
        # EVOLUTIONARY VERSION
        (sd,_) = assemblyevo(initpoolnet, rates0, maxits, cn,cm,ce,cpred, diverse, restrict_colonization);
        
        trophLevels = sd.trophLevels
        heatMapMax[primInd,secInd,repetition] = mean(maximum.(trophLevels))
        heatMapMean[primInd,secInd,repetition] = mean(mean.(trophLevels))

        jldsave("data/$(simulationName)/$(paramName)=($(rPrimExt),$(rSecExt))_repet=$repetition.jld2",compress; simulationData = sd, rates0)
    end

    jldsave("data/$(simulationName)/results.jld2",compress; heatMapMax, heatMapMean)
    
    plotlyjs()
    maxPlot = Plots.heatmap(paramVals, paramVals, dropdims(mean(heatMapMax,dims=3),dims=3),
        size = (1280,720), xlabel = "primary extinction rate", ylabel = "secondary extinction rate",
        title = "Average maximal trophic level in itterations 9500 to 10000")
    meanPlot = Plots.heatmap(x = paramVals, y = paramVals, z=dropdims(mean(heatMapMean,dims=3),dims=3),
        size = (640,480), xlabel = "primary extinction rate", ylabel = "secondary extinction rate",
        title = "Average mean trophic level in itterations 9500 to 10000")

    display(maxPlot)
    display(meanPlot)

    Plots.savefig(maxPlot,"data/$simulationName/plots/maxTrophLevelPlot.png");
    Plots.savefig(meanPlot,"data/$simulationName/plots/meanTrophLevelPlot.png");
end

simulation()