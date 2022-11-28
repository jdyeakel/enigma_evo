if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/enigma_evo/EnigmaEvo/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
end

#load data
filename = "data/nicheconstruction/simsettings.jld";
namespace_settings = smartpath(filename);
@load namespace_settings S lambda SSprobs SOprobs OOprobs e_t n_t maxits probmut cn ce cpred cmvec lcmvec reps parametervec its;

maxsprich = SharedArray{Int64}(lcmvec,reps);
maxsprichtime = SharedArray{Float64}(lcmvec,reps);
collapsetime = SharedArray{Float64}(lcmvec,reps);

@time @sync @distributed for i=1:its
    
    cm_pos = parametervec[i,1];
    cm = cmvec[cm_pos];

    r = parametervec[i,2];

    #load data
    filename = "data/nicheconstruction/simdata.jld";
    indices = [cm_pos,r];
    namespace = smartpath(filename,indices);
    
    # @save namespace mass1 mass2 epsilonvec clock popstate toothdrop state toothlength1 toothlength2;
    
    @load namespace intm eb nb nb0 mb SSpwp SOpwp sprich rich mstrength evolvedstrength clock CID intm_evo mutstep freqe freqn events eb_evo nb_evo nb0_evo mb_evo;

    #Calculate max sprich
    maxsprich[cm_pos,r] = findmax(sprich)[1];
    maxsprichtime[cm_pos,r] = clock[findmax(sprich)[2]];

    #Calculate collapse time
    # minsprich = findall(!iszero,reverse(diff(sprich)))[1]
    # collapsetime[cm_pos,r] = clock[findall(x->x==0,sprich[10:maxits])[1]];
    collapsetime[cm_pos,r] = clock[maxits - findall(!iszero,reverse(diff(sprich)))[1]]

end

decaytime = collapsetime .- maxsprichtime;

avgmaxsprich = vec(mean(maxsprich,dims=2));
avgmaxsprichtime = vec(mean(maxsprichtime,dims=2));
avgcollapsetime = vec(mean(collapsetime,dims=2));
avgdecaytime = vec(mean(decaytime,dims=2));

scatterplot(avgmaxsprich)
scatterplot(avgmaxsprichtime)
scatterplot(avgcollapsetime)
scatterplot(avgdecaytime)

filename = "figures/fig_nicheconstruction_richness.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
# pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pal = brewer.pal(3,'Set1')
pdf($namespace,width=6,height=5)
plot($cmvec,$avgmaxsprich,xlab='Niche construction strength',ylab='Peak richness',pch=16)
dev.off()
"""


filename = "figures/fig_nicheconstruction_collapse.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
# pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pal = brewer.pal(3,'Set1')
pdf($namespace,width=6,height=5)
plot($cmvec,$avgcollapsetime,xlab='Niche construction strength',ylab='Time to collapse',pch=16)
dev.off()
"""

