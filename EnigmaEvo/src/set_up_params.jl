S::Int64 = 200;
maxits::Int64 = 1000;
const SOprobs = (   #as link types are mutually exclusive pn + pe + pm(function of lambda) <= 1 must be fulfilled!!!
p_n=0.002, #0.002,
p_e=0.01 #0.01
);
const SSmult = 1.0; OOmult = 0.0;
const SSprobs = (p_n = SSmult .* SOprobs.p_n , p_e = SSmult .* SOprobs.p_e);
const OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);

#Competitive gain of a make
const cm = Float64(pi);
#Competitive gain of a need
const cn = exp(1);
#Competitive loss of an eat
const ce = sqrt(2);
#Competitive loss from a predator
const cpred = 1.;   #everywhere else cf but not here for compatibility with older version

#expected objects per species
const lambda = 0.1;
const e_t = 0.; #always set to 0

const n_t = 1.; #always set to 1

#rc = Colonization rates
#rprimext = Local primary species extinction rate
#rsecext = Local secondary species extinction rate
#reo = Local object extinction rate
#revo = Evolutionary rate
#rext = Global extinction rate
rates0 = (rc = 1., rprimext = 1., rsecext = 4., reo = 1., revo = 0.05, rext = 0.0);

#Turn diversification dynamic on or off
# 0 = off
# 1 = on
diverse::Int = 0;

logging::Bool = true;