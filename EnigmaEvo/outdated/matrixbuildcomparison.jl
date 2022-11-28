if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
end

# loadfunc = include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl");
# @everywhere include("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/assembly2.jl")

probmut = 0.5;

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
lambda = 0.0;
athresh = 0; #always set to 0
nthresh = 1.0; #always set to 1
MaxN = convert(Int64,floor(S + S*lambda));

#link density
ld = Array{Float64}(undef,1000,2);
tm = Array{Float64}(undef,1000,2);

for i=1:1000


    edgelist,sID,oID = intmatrixv5(S,lambda,SSprobs,SOprobs,OOprobs);
    intm, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);
    e_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(intm);

    

    #Calculate stats
    #1 Link density
    ld[i,1] = sum(e_b)/S;
    ld[i,2] = length(findall(x->x==1,edgelist[:,3]))/S;
end

