simulation_name = "vary_cn_no_engineers_20_000_it_restrict_col";

avg_freqe, avg_freqe_pool, avg_freqn, avg_freqn_pool, avg_sprich, avg_pool  = 
    average_time_series(["freqe","freqe_pool","freqn","freqn_pool","sprich","pool"],
                        simulation_name,"cn",0:0.5:5,50)

jldsave("data/$(simulation_name)/averages.jld2",true;avg_freqe, avg_freqe_pool, avg_freqn, avg_freqn_pool, avg_sprich, avg_pool)

avg_freqe, avg_freqe_pool, avg_freqn, avg_freqn_pool, avg_sprich, avg_pool = 
    load("EnigmaEvo/data/$(simulation_name)/averages.jld2", "avg_freqe", "avg_freqe_pool",
        "avg_freqn", "avg_freqn_pool", "avg_sprich", "avg_pool");

itterations = 1:20_000

lineplot(itterations,hcat(avg_sprich[:,[1,end]],avg_pool[:,[1,end]]),xlabel="itteration",ylabel="species richness",title="Species richness of new model")
lineplot(itterations,hcat(avg_freqe[:,[1,end]],avg_freqe_pool[:,[1,end]]),xlabel="itteration",ylabel="freq eat",title="Frequency of eat interactions of new model")
lineplot(itterations,hcat(avg_freqn[:,[1,end]],avg_freqn_pool[:,[1,end]]),xlabel="itteration",ylabel="freq need",title="Frequency of need interactions of new model")

freqe, freqn,freqe_pool,freqn_pool, sprich,pool = load("EnigmaEvo/data/vary_cn_no_engineers/cn=5.0_repet=2.jld2", "freqe", "freqn","freqe_pool","freqn_pool", "sprich","pool");
itterations = 1:20000

lineplot(itterations,hcat(sprich,pool),xlabel="itteration",ylabel="species richness",title="Species richness of new model")
lineplot(itterations,hcat(freqe,freqe_pool),xlabel="itteration",ylabel="freq eat",title="Frequency of eat interactions of new model")
lineplot(itterations,hcat(freqn,freqn_pool),xlabel="itteration",ylabel="freq need",title="Frequency of need interactions of new model")


paramName = "nBasalRes"
simulationName = "vary_$(paramName)_50_000_it";        #specify the name of the simulation

specRich_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "species richness",legend = :topright);
meanEats_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "average amount of eat interactions",legend = :topright);
meanNeeds_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "average amount of need interactions",legend = :topright);

finalSpecRichPlt = plot(xlabel="number of basal resources",ylabel="mean spec richness in last 500 itterations", legend=nothing);

param_vals = 11:10:101
for nBasalRes in param_vals
    for rep in 1:10
        filename = "$(paramName)=$(nBasalRes)_repet=$(rep).jld2"
        sd = load("data/$(simulationName)/$(filename)", "simulationData");
        if  sd.clock[end] != 0. && maximum(sd.clock) < 1000
            plot!(specRich_plt,sd.clock,sd.specRich,color=RGB(nBasalRes/param_vals[end],0,0),label = "$(paramName) = $(nBasalRes)");
            plot!(meanEats_plt,sd.clock,sd.meanEats,color=RGB(nBasalRes/param_vals[end],0,0),label = "$(paramName) = $(nBasalRes)");
            plot!(meanNeeds_plt,sd.clock,sd.meanNeeds,color=RGB(nBasalRes/param_vals[end],0,0),label = "$(paramName) = $(nBasalRes)");
            scatter!(finalSpecRichPlt,[nBasalRes], [mean(sd.specRich[45_000:end])]);
        end
    end
end

Plots.savefig(specRich_plt,"data/$simulationName/plots/specRich_plt.png");
Plots.savefig(meanEats_plt,"data/$simulationName/plots/meanEats_plt.png");
Plots.savefig(meanNeeds_plt,"data/$simulationName/plots/meanNeeds_plt.png");
Plots.savefig(finalSpecRichPlt,"data/$simulationName/plots/finalSpecRichPlt.png");

display(specRich_plt)
display(meanNeeds_plt)
display(meanEats_plt)

display(finalSpecRichPlt)



node = getnode(phyloTree,getparent(phyloTree,"51v2"))
node.data["evolution"]


maxTrophLevels = maximum.(sd.trophLevels)
findmax(maximum.(sd.trophLevels))

col_2560 = extraData.oldColonies[2560]
pool_2560 = extraData.oldPools[2560]

colnet_1950 = recreatecolnetdiverse(sd,1950)
colnet_1950.spec

eatMatrix = ENIgMaGraphs.convertToEatMatrixNonReduced(col_2560)
eatMatrixReduced = ENIgMaGraphs.convertToEatMatrix(col_3770)
R"""
    library(MASS)  
    library(NetIndices)
    rtl<-TrophInd(t($eatMatrix))
"""
@rget rtl;
trophLevels = rtl[:,:TL]
findmax(trophLevels)

"$(sort(trophLevels))"

inds = sort(collect(col_3770.spec))


eatMatrix[inds,inds] == eatMatrixReduced

function pathes_to_basal(net,id)
    pathLengths = Dict{Int,Vector{Int}}()
    function followPath(path)
        currId = path[end]
        if haskey(pathLengths,currId)
            push!(pathLengths[currId], length(path))
        else
            pathLengths[currId] = [length(path)]
        end

        for specId in net[currId].feed
            if !(specId in path)
                followPath(push!(path, specId))
            end
        end
    end

    for basalResId in net.basalRes
        for specId in net[basalResId].feed
            followPath([specId])
        end
    end
    ret = []
    for path in values(pathLengths)
        if path[end] == id
            push!(ret,path)
        end
    end
    return ret
end

pathes = pathes_to_basal(col_2560,208)

ENIgMaGraphs.getConnectedSpec(col_2560)
col_2560.spec
sd.nPrimExtSpec[2560]

mean(maximum.(sd.trophLevels)[900:end])