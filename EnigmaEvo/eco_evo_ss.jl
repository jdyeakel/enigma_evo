if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
end

# loadfunc = include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
# @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/assembly2.jl")

S = 200;
maxits = 4000;
SOprobs = (
p_n=0.002,
p_a=0.01
);
SSmult = 1.0; OOmult = 0.0;
SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);

cn = sqrt(2);
ce = pi;
cp = 1;


#expected objects per species
lambda = 0;
athresh = 0;
nthresh = 1.0;
MaxN = convert(Int64,floor(S + S*lambda));

probmutvec = collect(0:0.01:0.9)

@sync @distributed for i=1:length(probmutvec)

    probmut = probmutvec[i];


    intm, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);
    intm_origin = copy(intm);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(intm);

    sprich,rich,clock,CID,intm_evo,mutstep,freqe,freqn,events = assemblyevo(
        intm,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,lambda,
        athresh,nthresh,maxits,probmut,cn,ce,cp);

    a_b_origin,
    n_b_origin,
    i_b_origin,
    m_b_origin,
    n_b0_origin,
    sp_v,
    int_id = preamble_defs(intm_origin);
    a_b_evo,
    n_b_evo,
    i_b_evo,
    m_b_evo,
    n_b0_evo,
    sp_v,
    int_id = preamble_defs(intm_evo);

    #Rerun *evolved* system with eco-assembly
    sprich_evoeco,rich_evoeco,clock_evoeco,CID_evoeco,events_evoeco = assemblyeco(
    intm_evo,a_b_evo,n_b_evo,i_b_evo,m_b_evo,n_b0_evo,sp_v,int_id,lambda,athresh,nthresh,maxits,cn,ce,cp);

    #Rerun from initial state with eco-assembly
    sprich_eco,rich_eco,clock_eco,CID_eco,events_eco = assemblyeco(
    intm_origin,a_b_origin,n_b_origin,i_b_origin,m_b_origin,n_b0_origin,sp_v,int_id,lambda,athresh,nthresh,maxits,cn,ce,cp);
    
    #Save data
    if i == 1
        #Save data
        filename = "data/ss/sim_settings.jld";
        namespace = smartpath(filename);
        @save namespace lambda cn ce cp maxits probmutvec;
    end
    
    filename = "data/ss/cid.jld";
    indices = [i];
    namespace = smartpath(filename,indices);
    @save namespace intm_origin intm_evo sprich rich clock CID mutstep freqe freqn events sprich_evoeco rich_evoeco clock_evoeco CID_evoeco events_evoeco sprich_eco rich_eco clock_eco CID_eco events_eco;
    
end

