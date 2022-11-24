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
simulationName = "vary_$(paramName)";        #specify the name of the simulation

specRich_plt = plot(size = (4000,4000),xlabel = "clock time", ylabel = "species richness",legend = :topleft);
meanEats_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "average amount of eat interactions",legend = :topleft);
meanNeeds_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "average amount of need interactions",legend = :topleft);

param_vals = 1:10:101
for param in param_vals
    for rep in 1:2:9
        filename = "$(paramName)=$(param)_repet=$(rep).jld2"
        sd = load("data/$(simulationName)/$(filename)", "simulationData");
        plot!(sd.specRich_plt,sd.clock,sd.specRich,color=RGB(param/param_vals[end],0,0),label = "$(paramName) = $(param)");
        plot!(sd.meanEats_plt,sd.clock,sd.meanEats,color=RGB(param/param_vals[end],0,0),label = "$(paramName) = $(param)");
        plot!(sd.meanNeeds_plt,sd.clock,sd.meanNeeds,color=RGB(param/param_vals[end],0,0),label = "$(paramName) = $(param)");
        end
end

Plots.savefig(specRich_plt,"data/$simulationName/specRich_plt.png");
Plots.savefig(meanEats_plt,"data/$simulationName/meanEats_plt.png");
Plots.savefig(meanNeeds_plt,"data/$simulationName/meanNeeds_plt.png");

display(specRich_plt)
display(meanNeeds_plt)
display(meanEats_plt)



node = getnode(phyloTree,getparent(phyloTree,"51v2"))
node.data["evolution"]