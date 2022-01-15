function strengthcalc(nb0mut::BitArray{2},ebmut::BitArray{2},cid::Array{Int64},spmut::Int64,cn::Float64,ce::Float64,cp::Float64)::Float64
    spstrength = (cn*sum(nb0mut[spmut,cid])) .- (ce*sum(ebmut[spmut,:])) .- (cp*sum(ebmut[cid,spmut]));
    return spstrength
end
