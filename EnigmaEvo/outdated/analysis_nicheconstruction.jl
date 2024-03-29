if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/enigma_evo/EnigmaEvo/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
end

S = 200;
maxits = 5000;
SOprobs = (
p_n=0.002, #0.002,
p_e=0.01 #0.01
);
SSmult = 1.0; OOmult = 0.0;
SSprobs = (p_n = SSmult .* SOprobs.p_n , p_e = SSmult .* SOprobs.p_e);
OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);

#Competitive gain of a make
# cm = Float64(pi);
#Competitive gain of a need
cn = exp(1);
#Competitive loss of an eat
ce = sqrt(2);
#Competitive loss from a predator
cpred = 1.;


#expected objects per species
lambda = 0.1;
e_t = 0.; #always set to 0
n_t = 1.; #always set to 1
# MaxN = convert(Int64,floor(S + S*lambda));

#rc = Colonization rates
#re = Local species extinction rate
#reo = Local object extinction rate
#revo = Evolutionary rate
#rext = Global extinction rate
rates0 = (rc = 1., re = 1., reo = 1., revo = 1., rext = 0.01);
#Turn diversification dynamic on or off
# 0 = off
# 1 = on
diverse = 0;




reps = 1000;
cmvec = collect(0.:Float64(pi)/20:5*Float64(pi));
lcmvec = length(cmvec);

parametervec = [repeat(collect(1:lcmvec),inner=reps) repeat(collect(1:reps),outer=lcmvec)];
its = size(parametervec)[1];

#save data
filename = "data/nicheconstruction/simsettings.jld";
namespace_settings = smartpath(filename);
@save namespace_settings S lambda SSprobs SOprobs OOprobs e_t n_t maxits probmut cn ce cpred cmvec lcmvec reps parametervec its diverse rates0;

@time @sync @distributed for i=1:its
    
    cm_pos = parametervec[i,1];
    cm = cmvec[cm_pos];

    r = parametervec[i,2];


    intm,eb,nb,nb0,mb,SSpwp,SOpwp = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

    # edgelist_origin,sID,oID = intmatrixv5(S,lambda,SSprobs,SOprobs,OOprobs);
    # length(findall(x->x==1,edgelist_origin[:,3]))/S^2
    # length(findall(x->x==2,edgelist_origin[:,3]))/S^2

    @time sprich,rich,pool,mstrength,evolvedstrength,clock,CID,intm_evo,mutstep,freqe,freqn,events = assemblyevo(rates0,S,intm,eb,nb,nb0,mb,e_t,n_t,maxits,cm,cn,ce,cpred,diverse); eb_evo,nb_evo,nb0_evo,mb_evo = intbool(intm_evo);


    #save data
    filename = "data/nicheconstruction/simdata.jld";
    indices = [cm_pos,r];
    namespace = smartpath(filename,indices);
    
    # @save namespace mass1 mass2 epsilonvec clock popstate toothdrop state toothlength1 toothlength2;
    
    @save namespace intm eb nb nb0 mb SSpwp SOpwp sprich rich mstrength evolvedstrength clock CID intm_evo mutstep freqe freqn events eb_evo nb_evo nb0_evo mb_evo;

end
