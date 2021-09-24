function potcol2(edgelist::Array{Int64},splist::Array{Int64},cid::Array{Int64},e_t::Float64,n_t::Float64)::Array{Int64}
    #Species not in the community
    notcid = setdiff(splist,cid);
    
    #COUNT POTENTIAL COLONIZERS
    trophiclinked_bool = Array{Bool}(undef,length(notcid));
    servicelinked_bool = Array{Bool}(undef,length(notcid));
    for i=1:length(notcid)
        interactor = notcid[i];
        #Is there AT LEAST one trophic interaction available in community?
        trophic_interactions= intersect(intfind(edgelist,interactor,1),cid);
        #NOTE: eat threshold could be built in here
        trophiclinked_bool[i] = length(trophic_interactions) > 0;

        #Are all need interactions serviced in the community?
        need_interactions = intfind(edgelist,interactor,2);
        #NOTE: n_threshold could be built in here
        needsnotfilled = sum(indexin(need_interactions,cid) .=== nothing);
        servicelinked_bool[i] = needsnotfilled == 0;
    end
    trophiclinked = notcid[trophiclinked_bool];
    servicelinked = notcid[servicelinked_bool];
    col = intersect(trophiclinked,servicelinked);
    return col
end
