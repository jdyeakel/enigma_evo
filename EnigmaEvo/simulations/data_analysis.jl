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
