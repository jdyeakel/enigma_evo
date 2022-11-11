mutable struct PhyloNode
    const id::Int
    const timestamp::Float64
    descendants::Union{(PhyloNode,PhyloNode),Nothing}
    PhyloNode(id,timestamp) = new(id,timestamp,nothing)
end

