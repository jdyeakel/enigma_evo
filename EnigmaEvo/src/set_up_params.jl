S::Int64 = 200;
maxits::Int64 = 4_000;

nBasalResources::Int = 15
const SOprobs = (   #as link types are mutually exclusive pn + pe + pm(function of lambda) <= 1 must be fulfilled!!!
p_n=0.002, #0.002,
p_e=0.01 #0.01
);
const SSmult = 1.0; OOmult = 0.0;
const SSprobs = (p_n = SSmult .* SOprobs.p_n , p_e = SSmult .* SOprobs.p_e);
const OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);

#Competitive gain of a make
cm::Float64 = Float64(pi);
#Competitive gain of a need
cn::Float64 = exp(1);
#Competitive loss of an eat
ce::Float64 = sqrt(2);
#Competitive loss from a predator
cpred::Float64 = 1.;   #everywhere else cf but not here for compatibility with older version

#expected objects per species
lambda::Float64 = 0.0;
e_t::Float64 = 0.; #always set to 0

n_t::Float64 = 1.; #always set to 1

#rc = Colonization rates
#rprimext = Local primary species extinction rate
#rsecext = Local secondary species extinction rate
#reo = Local object extinction rate
#revo = Evolutionary rate
#rext = Global extinction rate
rates0::NamedTuple{(:rc, :rprimext, :rsecext, :reo, :revo, :rext), NTuple{6, Float64}} =
    (rc = 1., rprimext = 10., rsecext = 30., reo = 1., revo = 0.035, rext = .005);

#Turn diversification dynamic on or off
# 0 = off
# 1 = on
diverse::Int = 1;

#if true only species that wouldnt go to secondary (structural) extionction immediatelly can colonize
#otherwise all species can colonize
restrict_colonization::Bool = true;

logging::Bool = true;
