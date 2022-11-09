function plot_simulation(simulation_data;offset=0,show=true)
    poolnet,colnet,sprich,rich,pool,mstrength,evolvedstrength,clock,CID,maxids,glob_ext_spec,mutstep,freqe,freqn,freqe_pool,freqn_pool,events = simulation_data

    sprich_plt = plot(clock, sprich, xlabel="clock time", ylabel="species richness", legend=false);
    freqe_plt = plot(clock, freqe, xlabel="clock time", ylabel="average amount of eat interactions", legend=false);
    freqn_plt = plot(clock, freqn, lxlabel="clock time", ylabel="average amount of need interactions", legend=false);

    ext_len_dist = get_extinction_size_distrib(sprich[(1+offset):end])
    ext_size_plt = bar(pairs(ext_len_dist), xlabel="length of extinction cascade", yaxis=:log, ylabel="log(probability)", legend=false);

    l = @layout [
        [grid(3,1)] d{0.5w}
        ]
    summary_plt = plot(sprich_plt,freqe_plt,freqn_plt,ext_size_plt, layout=l,legend=false);
    if show
        display(summary_plt)
    end
    return summary_plt
end