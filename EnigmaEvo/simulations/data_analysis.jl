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


param_name = "cn"
simulation_name = "vary_$(param_name)_no_engineers_rprimext=10";        #specify the name of the simulation

sprich_plt = plot(size = (4000,4000),xlabel = "clock time", ylabel = "species richness",legend = :topleft);
freqe_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "average amount of eat interactions",legend = :topleft);
freqn_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "average amount of need interactions",legend = :topleft);

param_vals = 0:2:20
for param in param_vals
    for rep in 1:2
        filename = "$(param_name)=$(param)_repet=$(rep).jld2"
        clock, freqe, freqn, sprich =
            load("data/$(simulation_name)/$(filename)", "clock","freqe", "freqn", "sprich");
        plot!(sprich_plt,clock,sprich,color=RGB(param/param_vals[end],0,0),label = "cn = $(param)");
        plot!(freqe_plt,clock,freqe,color=RGB(param/param_vals[end],0,0),label = "cn = $(param)");
        plot!(freqn_plt,clock,freqn,color=RGB(param/param_vals[end],0,0),label = "cn = $(param)");
        end
end

Plots.savefig(sprich_plt,"data/$simulation_name/sprich_plt.png");
Plots.savefig(freqe_plt,"data/$simulation_name/freqe_plt.png");
Plots.savefig(freqn_plt,"data/$simulation_name/freqn_plt.png");

display(sprich_plt)
display(freqn_plt)
display(freqe_plt)



node = getnode(phyloTree,getparent(phyloTree,"51v2"))
node.data["evolution"]