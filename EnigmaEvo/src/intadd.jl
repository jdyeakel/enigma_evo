function intadd(edgelist::Array{Int64},spp::Array{Int64},edgetype)::Array{Int64}
    #Input: edge list, interactor, interactee
    #Output: updated edgelist (if interaction does not already exist)

    # #Is interactor on list? (FOR ANY TYPE OF INTERACTION)
    #NOTE: This will only check if an indegree exists, and not outdegree - keep in mind
    pos_interactor = edgelist[:,1] .== spp[1];
    # #Is interactee on list?
    pos_interactee = edgelist[:,2] .== spp[2];

    #Does this interaction exist?
    test = any(pos_interactor .* pos_interactee);
    #If not, then add the interaction
    if test == false
        vecset = [spp[1] spp[2]; spp[2] spp[1]];
        newset = [vecset[1:length(edgetype),:] edgetype']

        edgelist_update = vcat(edgelist,newset);
        return edgelist_update
    else
        return edgelist
    end
    
end
