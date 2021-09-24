function potobextinct2(edgelist::Array{Int64},ocid::Array{Int64},cid::Array{Int64})::Array{Int64}
     #COUNT POTENTIAL EXTINCT SPECIES
    makerlinked = falses(length(ocid));
    for i=1:length(ocid)
        object = ocid[i];
        makers = length(intersect(intfind_out(edgelist,object,3),cid));
        makerlinked[i] = makers == 0;
    end
    obext = ocid[makerlinked];
    return obext
end

