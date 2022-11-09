function potsecextinct2(edgelist::Array{Int64},spcid::Array{Int64},cid::Array{Int64},cn::Float64,ce::Float64,cp::Float64)::Array{Int64}
     #COUNT POTENTIAL EXTINCT SPECIES
    
    #1) Calculate strength of each species (not objects)
    spstrength = Array{Float64}(undef,length(spcid));
    for i = 1:length(spcid)
        interactor = spcid[i];
        spstrength[i] = strength(edgelist,interactor,cid,cn,ce,cp);
        # numberofneeds = length(intfind(edgelist,interactor,2)); #TOTAL Needs
        # numberofeats = length(intersect(intfind(edgelist,interactor,1),cid)); #REALIZED Eats
        # numberofpreds = length(intersect(intfind_out(edgelist,interactor,1),cid));
        # spstrength[i] = cn*numberofneeds - ce*numberofeats - cp*numberofpreds;
    end
    
    #2) Is the strength of any species less than the strength of a complete competitor?
    # DONT COUNT PRIMARY RESOURCE
    # prext_comp = Array{Bool}(undef,length(spcid)) .* false;
    prext_comp = falses(length(spcid));
    for i=1:length(spcid)
        interactor = spcid[i];
        #Trophic interactions with trophic level ABOVE the primary resource
        trophic_interactions = setdiff(intersect(intfind(edgelist,interactor,1),cid),1);
        for j=1:length(trophic_interactions)
            potential_competitors = setdiff(intersect(intfind_out(edgelist,trophic_interactions[j],1),cid),interactor);
            for k = 1:length(potential_competitors)
                #Position of competitor in spcid
                cpos = indexin(potential_competitors[k],spcid)[1];
                #List of competitors realized foods
                c2 = setdiff(intersect(intfind(edgelist,potential_competitors[k],1),cid),1);
                #If list of competitors realized foods is same AND its strength is greater, then interactor subject to extinction
                prext_comp[i] = (trophic_interactions == c2) & (spstrength[i] < spstrength[cpos]);
            end
        end
    end
    spext2 = spcid[prext_comp];
    return spext2
end
