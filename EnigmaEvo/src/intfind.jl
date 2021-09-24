function intfind(edgelist::Array{Int64},interactor::Int64,edgetype::Int64)::Array{Int64}
    #Subset of edgelist for interactors

    if edgetype > 0
        interactees = edgelist[((edgelist[:,1] .== interactor) .& (edgelist[:,3] .== edgetype)),2];
    else
        interactees = edgelist[(edgelist[:,1] .== interactor),2];
    end


    # interactorpos = findall(x->x==interactor,edgelist[:,1]);
    # intedges = edgelist[interactorpos,:];

    # interactees = intedges[findall(x->x==edgetype,intedges[:,3]),2]

    splist = unique(interactees);

    return sort(splist)
end
