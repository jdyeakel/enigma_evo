Base.eltype(::Type{<:AbstractENIgMaGraph}) = Int

Graphs.is_directed(::Type{<:AbstractENIgMaGraph}) = true
Graphs.is_directed(::AbstractENIgMaGraph) = true

@enum InteractionType begin
    eatInteraction=1
    needInteraction=2
    makeInteraction=3
end

struct ENIgMaEdge{intType}
    src::Int
    dst::Int
end

struct InteractionGraph{intType} <: AbstractENIgMaGraph
    g::ENIgMaGraph
end

EatGraph = InteractionGraph{eatInteraction}
NeedGraph = InteractionGraph{needInteraction}
MakeGraph = InteractionGraph{makeInteraction}


Graphs.edges(g::EatGraph) = 
    (ENIgMaEdge{eatInteraction}(idSrc, idDest) 
        for idSrc in g.g.spec for idDest in g.g[idSrc].eat)
    
Graphs.edges(g::NeedGraph) = 
    (ENIgMaEdge{needInteraction}(idSrc, idDest) 
        for idSrc in g.g.spec for idDest in g.g[idSrc].need )

Graphs.edges(g::MakeGraph) = 
    (ENIgMaEdge{makeInteraction}(idSrc, idDest) 
        for idSrc in g.g.spec  for idDest in g.g[idSrc].make )
            
asEatGraph(g::ENIgMaGraph) = EatGraph(g)
asNeedGraph(g::ENIgMaGraph) = NeedGraph(g)
asMakeGraph(g::ENIgMaGraph) = MakeGraph(g)
asInteractionGraph(intType, g::ENIgMaGraph) = InteractionGraph{intType}(g) 

Graphs.edges(g::ENIgMaGraph) =
    Iterators.flatten((edges(asEatGraph(g)), edges(asNeedGraph(g)), edges(asMakeGraph(g))))

Graphs.edgetype(g::AbstractENIgMaGraph) = ENIgMaEdge

Graphs.edgetype(g::InteractionGraph{intType}) where intType = ENIgMaEdge{intType}

function Graphs.has_edge(g::ENIgMaGraph, srcId, destId)
    srcId in g || return false
    srcNode = g[srcId]
    destId in srcNode.eat || destId in srcNode.need || destId in srcNode.make
end 

function Graphs.has_edge(g::EatGraph,srcId, destId)
    srcId in g.g || return false
    destId in g.g[srcId].eat
end

function Graphs.has_edge(g::NeedGraph,srcId, destId)
    srcId in g.g || return false
    destId in g.g[srcId].need
end

function Graphs.has_edge(g::MakeGraph,srcId, destId)
    srcId in g.g || return false
    destId in g.g[srcId].make
end

Graphs.has_vertex(g::ENIgMaGraph,id) = id in g
Graphs.has_vertex(g::InteractionGraph{intType},id) where intType = id in g.g

Graphs.inneighbors(g::EatGraph,id) = collect(g.g[id].feed)

function Graphs.inneighbors(g::NeedGraph,destId)
    destId in g.g || throw(KeyError(destId))
    [srcId for srcId in g.g.spec if destId in g.g[srcId].need]
end

Graphs.inneighbors(g::MakeGraph,id) = g.g.hasspec[id] ? [] : collect(g.g[id].need)    #only needs of modifiers are incoming makes

function Graphs.inneighbors(g::ENIgMaGraph,destId)
    destId in g || throw(KeyError(destId))
    destNode = g[destId]
    [
        collect(destNode.feed) 
        [srcId for srcId in g.spec if destId in g[srcId].need]
        (g.hasspec[destId] ? [] : collect(destNode.need))
    ]
end

Graphs.outneighbors(g::EatGraph,id) = collect(g.g[id].eat)

function Graphs.outneighbors(g::NeedGraph,id)
    if g.g.hasspec[id] #needs of modifiers arent proper needs
        return [destId for destId in g.g[id].need if destId in g.g] #only needs to vertices in network count
    else
        id in g.g || throw(KeyError(id))
        return []   
    end
end

Graphs.outneighbors(g::MakeGraph,id) = collect(g.g[id].make)

Graphs.outneighbors(g::ENIgMaGraph,id) =
    [
        outneighbors(asEatGraph(g),id)
        outneighbors(asNeedGraph(g),id)
        outneighbors(asMakeGraph(g),id)
    ]


Graphs.ne(g::EatGraph) = 
    sum( length(spec.eat) for spec in (g.g[specId] for specId in g.g.spec ))
Graphs.ne(g::NeedGraph) = 
    sum( sum(g.g.hasv[collect(spec.need)]) for spec in (g.g[specId] for specId in g.g.spec )) #only count needs actually realized in the network
Graphs.ne(g::MakeGraph) = 
    sum( length(spec.make) for spec in (g.g[specId] for specId in g.g.spec ))
Graphs.ne(g::ENIgMaGraph) = 
    sum( length(spec.eat) + sum(g.hasv[collect(spec.need)]) + length(spec.make)
        for spec in (g[specId] for specId in g.spec ))

Graphs.nv(g::InteractionGraph{intT}) where intT = length(g.g)
Graphs.nv(g::ENIgMaGraph) = length(g)

Graphs.vertices(g::InteractionGraph{intT}) where intT = keys(g.g.vert)
Graphs.vertices(g::ENIgMaGraph) = keys(g.vert)



