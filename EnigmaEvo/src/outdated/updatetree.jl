function updatetree(treenames,tree,newid,spmut)
    push!(treenames,newid);
    numsp = size(tree)[1];
    newtree = spzeros(numsp+1,numsp+1);
    newtree[1:numsp,1:numsp] = tree;
    newtree[(numsp+1),spmut] = 1;
    return(treenames,newtree)
end
