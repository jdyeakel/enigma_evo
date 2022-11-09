function potextinct2(edgelist::Array{Int64},spcid::Array{Int64},cid::Array{Int64})::Array{Int64}
     #COUNT POTENTIAL EXTINCT SPECIES
    #1) By not fulfilling Eat/Need thresholds

    #COUNT POTENTIAL EXTINCT
    trophicunlinked_bool = Array{Bool}(undef,length(spcid));
    serviceunlinked_bool = Array{Bool}(undef,length(spcid));
    for i=1:length(spcid)
        interactor = spcid[i];
        #Is there AT LEAST one trophic interaction available in community?
        trophic_interactions= intersect(intfind(edgelist,interactor,1),cid);
        #NOTE: eat threshold could be built in here
        trophicunlinked_bool[i] = length(trophic_interactions) == 0;

        #Are all need interactions serviced in the community?
        need_interactions = intfind(edgelist,interactor,2);
        #NOTE: n_threshold could be built in here
        needsnotfilled = sum(indexin(need_interactions,cid) .=== nothing);
        serviceunlinked_bool[i] = needsnotfilled != 0;
    end
    trophicunlinked = spcid[trophicunlinked_bool];
    serviceunlinked = spcid[serviceunlinked_bool];
    spext1 = union(trophicunlinked,serviceunlinked);

    return spext1
    
end
