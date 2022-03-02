if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/enigma_evo/EnigmaEvo/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
end


#Probability of event-mutation
probmut = 0.5;
#This is the threshold for diversification
div_t = 1.0

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
e_t = 0.; #always set to 0
n_t = 1.; #always set to 1
# MaxN = convert(Int64,floor(S + S*lambda));


reps = 100;
lambdavec = collect(0.:0.1:2.);
llambdavec = length(lambdavec);
cmvec = collect(0.:Float64(pi)/10:5*Float64(pi));
lcmvec = length(cmvec);

parametervec1 = [repeat(collect(1:llambdavec),inner=lcmvec) repeat(collect(1:lcmvec),outer=llambdavec)];
parametervec = [repeat(parametervec1,outer=reps) repeat(collect(1:reps),inner=llambdavec*lcmvec)];
its = size(parametervec)[1];

#save data
filename = "data/eng_nicheconstruction/simsettings.jld";
namespace_settings = smartpath(filename);
@save namespace_settings S lambdavec llambdavec SSprobs SOprobs OOprobs e_t n_t maxits probmut cn ce cpred cmvec lcmvec reps parametervec its;

@time @sync @distributed for i=1:its
    
    lambda_pos = parametervec[i,1];
    lambda = lambdavec[lambda_pos];

    cm_pos = parametervec[i,2];
    cm = cmvec[cm_pos];

    r = parametervec[i,3];


    intm,eb,nb,nb0,mb,SSpwp,SOpwp = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

    # edgelist_origin,sID,oID = intmatrixv5(S,lambda,SSprobs,SOprobs,OOprobs);
    # length(findall(x->x==1,edgelist_origin[:,3]))/S^2
    # length(findall(x->x==2,edgelist_origin[:,3]))/S^2

    @time sprich,rich,mstrength,evolvedstrength,clock,CID,intm_evo,mutstep,freqe,freqn,events = assemblyevo(S,intm,eb,nb,nb0,mb,e_t,n_t,maxits,probmut,cm,cn,ce,cpred); eb_evo,nb_evo,nb0_evo,mb_evo = intbool(intm_evo);


    #save data
    filename = "data/eng_nicheconstruction/simdata.jld";
    indices = [cm_pos,r];
    namespace = smartpath(filename,indices);
    
    # @save namespace mass1 mass2 epsilonvec clock popstate toothdrop state toothlength1 toothlength2;
    
    @save namespace intm eb nb nb0 mb SSpwp SOpwp sprich rich mstrength evolvedstrength clock CID intm_evo mutstep freqe freqn events eb_evo nb_evo nb0_evo mb_evo;

end
