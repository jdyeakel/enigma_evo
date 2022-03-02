function strengthcalc(nb0::BitArray{2},eb::BitArray{2},mb::BitArray{2},cid::Array{Int64},sp::Int64,cm::Float64,cn::Float64,ce::Float64,cpred::Float64)::Float64
    spstrength = (cm*sum(mb[sp,cid]))  .+ (cn*sum(nb0[sp,cid])) .- (ce*sum(eb[sp,:])) .- (cpred*sum(eb[cid,sp]));
    return spstrength
end
