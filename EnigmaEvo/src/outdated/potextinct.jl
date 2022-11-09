function potextinct(spcid,cid,e_b,n_b0)
     #COUNT POTENTIAL EXTINCT SPECIES
    #1) By not fulfilling Eat/Need thresholds
    #Re-calculate athresh and nthresh (> athresh; >= nthresh)
    a_fill = ((sum(e_b[spcid,[1;cid]],dims=2)./sum(e_b[spcid,:],dims=2)) .> athresh)[:,1];
    prop_n = sum(n_b0[spcid,[1;cid]],dims=2)./sum(n_b0[spcid,:],dims=2);
    #If there are no 'need' interactions, proportion filled is always 1, and will always pass
    prop_n[isnan.(prop_n)] .= 1;
    n_fill = (prop_n .>= nthresh)[:,1];
    #Build a list of those that pass each individual test
    a_pass = spcid[a_fill];
    n_pass = spcid[n_fill];
    #which species pass?
    survivors = intersect(a_pass,n_pass);
    spext1 = setdiff(spcid,survivors);
    
    return spext1
end
