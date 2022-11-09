function secexteval(spcid,cid,eb,nb0,e_t,n_t)
    e_fill = ((sum(eb[spcid,[1;cid]],dims=2)./sum(eb[spcid,:],dims=2)) .> e_t)[:,1];
    prop_n = sum(nb0[spcid,[1;cid]],dims=2)./sum(nb0[spcid,:],dims=2);
    #If there are no 'need' interactions, proportion filled is always 1, and will always pass
    prop_n[isnan.(prop_n)] .= 1;
    n_fill = (prop_n .>= n_t)[:,1];
    #Build a list of those that pass each individual test
    e_pass = spcid[e_fill];
    n_pass = spcid[n_fill];
    #which species pass?
    survivors = intersect(e_pass,n_pass);
    secext = setdiff(spcid,survivors);
    return secext
end