if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/enigma_evo/EnigmaEvo/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
end

#load data
filename = "data/eng_nicheconstruction/simsettings.jld";
namespace_settings = smartpath(filename);
@load namespace_settings S lambdavec llambdavec SSprobs SOprobs OOprobs e_t n_t maxits probmut cn ce cpred cmvec lcmvec reps parametervec its diverse;

maxsprich = SharedArray{Int64}(llambdavec,lcmvec,reps);
maxsprichtime = SharedArray{Float64}(llambdavec,lcmvec,reps);
collapsetime = SharedArray{Float64}(llambdavec,lcmvec,reps);

@time @sync @distributed for i=1:its
    
    
    lambda_pos = parametervec[i,1];
    lambda = lambdavec[lambda_pos];

    cm_pos = parametervec[i,2];
    cm = cmvec[cm_pos];

    r = parametervec[i,3];

    #load data
    filename = "data/eng_nicheconstruction/simdata.jld";
    indices = [lambda_pos,cm_pos,r];
    namespace = smartpath(filename,indices);
    
    # @save namespace mass1 mass2 epsilonvec clock popstate toothdrop state toothlength1 toothlength2;
    
    @load namespace intm eb nb nb0 mb SSpwp SOpwp sprich rich mstrength evolvedstrength clock CID intm_evo mutstep freqe freqn events eb_evo nb_evo nb0_evo mb_evo;

    #Calculate max sprich
    maxsprich[lambda_pos,cm_pos,r] = findmax(sprich)[1];
    maxsprichtime[lambda_pos,cm_pos,r] = clock[findmax(sprich)[2]];

    #Calculate collapse time
    # minsprich = findall(!iszero,reverse(diff(sprich)))[1]
    # collapsetime[cm_pos,r] = clock[findall(x->x==0,sprich[10:maxits])[1]];

    if length(findall(!iszero,reverse(diff(sprich)))) > 0
        collapsetime[lambda_pos,cm_pos,r] = clock[maxits - findall(!iszero,reverse(diff(sprich)))[1]];
    else
        collapsetime[lambda_pos,cm_pos,r] = clock[maxits];
    end

end

decaytime = collapsetime .- maxsprichtime;

# avgmaxsprich = mean(maxsprich,dims=3)[:,:,1];
# avgmaxsprichtime = mean(maxsprichtime,dims=3)[:,:,1];
# avgcollapsetime = mean(collapsetime,dims=3)[:,:,1];
# avgdecaytime = mean(decaytime,dims=3)[:,:,1];

#There are sometimes NANs, so a pain!
for i=1:llambdavec
    for j=1:lcmvec
        avgmaxsprich[i,j] = mean(maxsprich[i,j,findall(!isinf,maxsprich[i,j,:])]);
        avgmaxsprichtime[i,j] = mean(maxsprichtime[i,j,findall(!isinf,maxsprichtime[i,j,:])]);
        avgcollapsetime[i,j] = mean(collapsetime[i,j,findall(!isinf,collapsetime[i,j,:])]);
        avgdecaytime[i,j] = mean(decaytime[i,j,findall(!isinf,decaytime[i,j,:])]);
    end
end






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


filename = "figures/fig_eng_nicheconstruction_collapse.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
# pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))(50)
pal = brewer.pal(4,'Set1')
pdf($namespace,width=5,height=5)
plot($cmvec,$(avgcollapsetime[1,:]),xlab='Niche construction strength',ylab='Time to collapse',pch=16,col=pal[1],ylim=c(20,70))
points($cmvec,$(avgcollapsetime[2,:]),pch=16,col=pal[2])
points($cmvec,$(avgcollapsetime[6,:]),pch=16,col=pal[3])
points($cmvec,$(avgcollapsetime[11,:]),pch=16,col=pal[4])
dev.off()
"""

