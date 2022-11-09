function strength(edgelist::Array{Int64},interactor::Int64,cid::Array{Int64},cn::Float64,ce::Float64,cp::Float64)::Float64
    numberofneeds = length(intfind(edgelist,interactor,2)); #TOTAL Needs
    numberofeats = length(intersect(intfind(edgelist,interactor,1),cid)); #REALIZED Eats
    numberofpreds = length(intersect(intfind_out(edgelist,interactor,1),cid));
    spstrength = cn*numberofneeds - ce*numberofeats - cp*numberofpreds;
    return spstrength
end
