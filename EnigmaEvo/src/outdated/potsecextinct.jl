function potsecextinct(spcid,cid,e_b,n_b0,cn,ce,cp,minstrength)
     #COUNT POTENTIAL EXTINCT SPECIES
    #1) By not fulfilling Eat/Need thresholds
    #Re-calculate athresh and nthresh (> athresh; >= nthresh)
    if length(spcid) > 0

        # Build the strength matrix at each community state
        # Strength values change over time so need to ve updated
        # Only record strength values of species (hence the [1:length(spcid)])
        #NOTE: Needs won't change; Eats is based on POTENTIAL niche; Vuln changes per timestep
        strength = vec(cn*sum(n_b0[spcid,cid],dims=2)) .- vec(ce*sum(e_b[spcid,:],dims=2)) .- (vec(cp*sum(e_b[spcid,cid],dims=1))[1:length(spcid)]);

        cmatrix = Array{Float64}(undef,length(spcid),length(cid));
        for i=1:length(spcid)
            cmatrix[i,:] .= e_b[spcid[i],cid] * strength[i];
        end

        # cmatrix = e_b[spcid,cid] .* reshape(repeat(strength,outer=length(cid)),length(spcid),length(cid));

        # cmatrix = (e_b[spcid,cid]' * (Matrix{Float64}(I,length(strength),length(strength)) .* strength))'

        #'zero' entrees need to be lower than any possible strength
        #So they will be effectively ignored
        #-sqrt(2)*S-S is the theoretical min strength
        cmatrix[cmatrix.==0] .= minstrength;

        #Maximum strength values for each resource utilized
        cmax = findmax(cmatrix,dims=1)[1];

        prext_comp = trues(length(spcid));
        for i=1:length(spcid)

            #Don't count sun
            #catalogue prey for all species/objects in the system (excluding sun)
            ieats = Array{Bool}(e_b[spcid[i],cid]);
            #If you have >= the max strength for any of those prey, you stay
            #This means that a pure primary producer is not evaluated
            prext_comp[i] = any(ieats)*(any(strength[i] .>= cmax[ieats])==false);

        end
        spext2 = spcid[prext_comp];
    else
        spext2 = Array{Int64}(undef,0);
    end

    return spext2
end
