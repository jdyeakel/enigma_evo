function plot_simulation(simulation_data;offset=0,show=true)
    poolnet, colnet, phyloTree, sprich,rich,pool,mstrength,evolvedstrength,clock,CID,maxids,glob_ext_spec,mutstep,freqe,freqn,freqe_pool,freqn_pool,events = simulation_data

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

function getAttributeVector(tree, attributeName, eltype)
    root = first(nodefilter(isroot, tree))
    ret = Vector{eltype}()
    function local!(tree,ret,node)
        isleaf(tree,node) && return
        push!(ret, node.data[attributeName])
        for child in getchildren(tree, node)
            local!(tree,ret,getnode(tree,child))
        end
    end
    local!(tree,ret,root)
    ret
end

function plotPhylogeny(phyloTree;sorted=true)
    sorted && sort!(phyloTree)
    
    evoTypes = ["i:ie","i:in","i:im","i:ei","i:en","i:em","i:ni","i:ne","i:nm","i:mi","i:me","i:mn",
                      "o:ie","o:in","o:ei","o:en","o:ni","o:ne"];
    evos = getAttributeVector(phyloTree,"evolution",String);
    nEvos = length(evos) - 1 #as long as artificial root is still neccessary
    evosDistribution = Dict{String,Float64}(evo => count((x)->x==evo,evos)/nEvos for evo in evoTypes)
    display(bar(pairs(evosDistribution),
        yaxis="probability",
        title = "Distribution of evolution types",
        size = (1000,400)
    ))

    plot(phyloTree, 
        series_annotations=text.(evos, 6, :center, :center, :black,),
        size=(1920,20_000),
        markersize=10,
        markerstrokewidth = .2,
        markercolor=:white,
        markershape=:rect
    )
    #plot(phyloTree, marker_group = "evolution", series_annotations = text.("test" for i in 1:nnodes(phyloTree)), markersize = 15, markershape = :square, size=(1920,20_000))
    #plt = plot(phyloTree, marker_group = "evolution", c= :Accent_8, markerstrokewidth = .2, legend=:topleft, markersize = 5, markershape = :square, size=(1920,20_000))
end