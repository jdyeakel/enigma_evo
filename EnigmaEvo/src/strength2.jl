function strengthcalc(nb0::BitArray{2},eb::BitArray{2},cid::Array{Int64},sp::Int64,cn::Float64,ce::Float64,cp::Float64)::Float64
    spstrength = (cm*sum(mb[sp,cid]))  .+ (cn*sum(nb0[sp,cid])) .- (ce*sum(eb[sp,:])) .- (cp*sum(eb[cid,sp]));
    return spstrength
end
