if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
end

#Load data
filename = "data/ss/sim_settings.jld";
namespace = smartpath(filename);
@load namespace lambda cn ce cp maxits probmutvec;

ssevo = Array{Float64}(undef,length(probmutvec));
sseco = Array{Float64}(undef,length(probmutvec));
ssevoeco = Array{Float64}(undef,length(probmutvec));
num_mut = Array{Float64}(undef,length(probmutvec));

for i=1:length(probmutvec)
    #Load data
    filename = "data/ss/cid.jld";
    indices = [i];
    namespace = smartpath(filename,indices);
    @load namespace intm_origin intm_evo sprich rich clock CID mutstep freqe freqn events sprich_evoeco rich_evoeco clock_evoeco CID_evoeco events_evoeco sprich_eco rich_eco clock_eco CID_eco events_eco;
    
    num_mut[i] = length(findall(!isodd,vec(intm_origin .== intm_evo))) / (size(intm_origin)[1]*(size(intm_origin)[2]-1));
    
    ssevo[i] = mean(sprich[2000:maxits]);
    sseco[i] = mean(sprich_eco[2000:maxits]);
    ssevoeco[i] = mean(sprich_evoeco[2000:maxits]);
end

filename = "/figures/steadystate_v2.pdf"
namespace = smartpath(filename);
R"""
pdf($namespace,width=5,height=4)
plot($num_mut, $ssevo / $sseco,pch=16,ylim=c(0,1.5),xlab='Proportion of evolved interactions',ylab='Steady state ratio')
points($num_mut, $ssevoeco / $sseco,pch=16,col='red')
lines($num_mut,rep(1,length($num_mut)),lty=3)
dev.off()
"""

filename = "/figures/steadystate2.pdf"
namespace = smartpath(filename);
R"""
pdf($namespace,width=5,height=4)
plot($num_mut, $ssevo,pch=16,ylim=c(1,200))
points($num_mut, $ssevoeco,pch=16,col='red')
points($num_mut, $sseco,pch=16,col='blue')
dev.off()
"""
