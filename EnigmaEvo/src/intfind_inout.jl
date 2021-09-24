function intfind_inout(edgelist::Array{Int64},interactor::Int64,interactee::Int64)::Tuple{Int64,Int64}
    #NOTE: This finds the 'out' interactions
    
    #Subset of edgelist for interactors

    # interactees = edgelist[((edgelist[:,1] .== interactor) .& (edgelist[:,3] .== edgetype)),2];
    position_test = ((edgelist[:,1] .== interactor) .& (edgelist[:,2] .== interactee))
    position_test_true = findall(isone,position_test);
    if length(position_test_true) > 0
        position_int = position_test_true[1];
        edge_type = edgelist[position_int,3];
    else 
        position_int = 0;
        edge_type = 0;
    end

    # interactorpos = findall(x->x==interactor,edgelist[:,1]);
    # intedges = edgelist[interactorpos,:];

    # interactees = intedges[findall(x->x==edgetype,intedges[:,3]),2]


    return position_int,edge_type
end
