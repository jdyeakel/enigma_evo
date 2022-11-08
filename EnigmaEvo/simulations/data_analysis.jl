avg_freqe, avg_freqn, avg_sprich = load("EnigmaEvo/data/vary_cn_no_engineers_2000_it/averages.jld2", "avg_freqe", "avg_freqn", "avg_sprich");

clock = 1:2000

lineplot(clock,avg_sprich[:,[1,end]],xlabel="clock time",ylabel="species richness",title="Species richness of new model")
lineplot(clock,avg_freqe[:,[1,end]],xlabel="clock time",ylabel="freq eat",title="Frequency of eat interactions of new model")
lineplot(clock,avg_freqn[:,[1,end]],xlabel="clock time",ylabel="freq need",title="Frequency of need interactions of new model")

freqe, freqn,freqe_pool,freqn_pool, sprich,pool = load("EnigmaEvo/data/vary_cn_no_engineers/cn=5.0_repet=2.jld2", "freqe", "freqn","freqe_pool","freqn_pool", "sprich","pool");
clock = 1:20000

lineplot(clock,hcat(sprich,pool),xlabel="clock time",ylabel="species richness",title="Species richness of new model")
lineplot(clock,hcat(freqe,freqe_pool),xlabel="clock time",ylabel="freq eat",title="Frequency of eat interactions of new model")
lineplot(clock,hcat(freqn,freqn_pool),xlabel="clock time",ylabel="freq need",title="Frequency of need interactions of new model")

sprich_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "species richness",legend = :topleft);
freqn_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "average amount of eat interactions",legend = :topleft);
freqe_plt = plot(size = (1920,1080),xlabel = "clock time", ylabel = "average amount of need interactions",legend = :topleft);

rep = 2
for cn in 0:0.5:5
    filename = "cn=$(cn)_repet=$(rep).jld2"
    clock,freqe, freqn, sprich = load("EnigmaEvo/data/vary_cn_no_engineers_20_000_it_restrict_col/$(filename)", "clock","freqe", "freqn", "sprich");
    plot!(sprich_plt,clock,sprich,color=RGB(cn/5,0,0),label = "cn = $(cn)")
    plot!(freqe_plt,clock,freqe,color=RGB(cn/5,0,0),label = "cn = $(cn)")
    plot!(freqn_plt,clock,freqn,color=RGB(cn/5,0,0),label = "cn = $(cn)")
end

display(sprich_plt)
display(freqn_plt)
display(freqe_plt)

