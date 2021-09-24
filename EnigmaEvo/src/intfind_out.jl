function intfind_out(edgelist::Array{Int64},interactee::Int64,edgetype::Int64)::Array{Int64}
    #NOTE: This finds the 'out' interactions
    
    #Subset of edgelist for interactors

    # interactees = edgelist[((edgelist[:,1] .== interactor) .& (edgelist[:,3] .== edgetype)),2];

    interactors = edgelist[((edgelist[:,2] .== interactee) .& (edgelist[:,3] .== edgetype)),1]

    # interactorpos = findall(x->x==interactor,edgelist[:,1]);
    # intedges = edgelist[interactorpos,:];

    # interactees = intedges[findall(x->x==edgetype,intedges[:,3]),2]

    splist = unique(interactors);

    return sort(splist)
end
