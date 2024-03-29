# defines a Graph type with multiple edge types stored in seperate adjacency lists (Sets actually)

#optimization ideas: include strength calculations in add and del calls

#performance questions: is there a good graph representation out there? Sets good choice in vert? dict good choice in graph

module ENIgMaGraphs

# exported functions variables are visible outside of the module if we include the module with using
export ENIgMaGraph, ENIgMaVert, IdManager       #might be put in begin block? kinda didnt work
export setUpPool, assemblyevo, ENIgMaSimulationData
export adde!, addn!, addf!, addm!, dele!, deln!, delm!, delf!
export getNumSpec,getNumBasalRes, getNumMods 
export addSpec!,replacespec!,addMod!, addBasalRes!, delv!
export getprimext, getsecext!, getpotcolonizers!
export colonize!,mutate!
export getNextId!
export converttoENIgMaGraph, converttointeractionmat, convertToEatMatrix
export getTrophicLevels, recreatecolnetdiverse, getConnectedSpec

export InteractionType, ignoreInteraction, eatInteraction, needInteraction, makeInteraction
export AbstractENIgMaEvent, ColonizationEvent, PrimaryExtinctionEvent, SecondaryExtinctionEvent
export ObjectExtinctionEvent, GlobalExtinctionEvent, MutationEvent, isMutationType
export getDegrees, getDegreeDist, modifiedSDSimilarities

using Distributed
@everywhere using Graphs
@everywhere using Phylo
@everywhere using Distributions
@everywhere using RCall
@everywhere using OffsetArrays

const enlargementfactor = 1.3; #controls how much buffer is added if estsize has to be increased

#is used to manage vertex ids, so far only tracks the highest id and assigns new ids
mutable struct IdManager
    maxid::Int

    IdManager(offset) = new(offset)
end;

setMaxId!(idmanager,maxid) = idmanager.maxid = maxid;

"""
    getnextid!(idmanager::IdManager)

    Returns the next available vertex id.
"""
function getNextId!(idmanager::IdManager) #optimization idea: have IdManager save ids of globally extinct species to reassign them?
    idmanager.maxid += 1;
    return idmanager.maxid;
end

"""
Defines a vertex in the ENIgMa-model.
"""
mutable struct ENIgMaVert
    eat::Set{Int}
    need::Set{Int}
    make::Set{Int}
    feed::Set{Int}

    strength::Float64

    #constructors 
    ENIgMaVert() = new(Set{Int}(),Set{Int}(),Set{Int}(),Set{Int}(),NaN) #creates an empty vertex
    ENIgMaVert(eat,need,make,feed) = new(eat,need,make,feed,NaN)    #creates a vertex with the given interactions
end;


#if changed: these have not been used all the time----------------------------

"""
    adde!(v::ENIgMaVert,ide)

    Adds 'ide' to the set of eat interactions of 'v'.
"""
adde!(v::ENIgMaVert,ide) = push!(v.eat,ide);

"""
    addn!(v::ENIgMaVert,idn)

    Adds 'idn' to the set of need interactions of 'v'.
"""
addn!(v::ENIgMaVert,idn) = push!(v.need,idn);
"""
    addm!(v::ENIgMaVert,idm)

    Adds 'idm' to the set of make interactions of 'v'.
"""
addm!(v::ENIgMaVert,idm) = push!(v.make,idm);
"""
    addf!(v::ENIgMaVert,idf)

    Adds 'idf' to the set of feed/predator interactions of 'v'.
"""
addf!(v::ENIgMaVert,idf) = push!(v.feed,idf);

function adde!(predator,predatorId,prey,preyId)   #creates an eat interaction together with the reciprocal feed interaction
    adde!(predator,preyId);
    addf!(prey,predatorId);
end

function addm!(spec,specid,mod,modid) #creates a make interaction together with the reciprocal need interaction.
    addm!(spec,modid);
    addn!(mod,specid);
end

#same as add functions just for deletion
"""
    dele!(v::ENIgMaVert,ide)

    Deletes 'ide' from the set of eat interactions of 'v'.
"""
dele!(v::ENIgMaVert,ide) = setdiff!(v.eat,ide);
"""
    deln!(v::ENIgMaVert,idn)

    Deletes 'idn' from the set of need interactions of 'v'.
"""
deln!(v::ENIgMaVert,idn) = setdiff!(v.need,idn);
"""
    delm!(v::ENIgMaVert,idm)

    Deletes 'idm' from the set of make interactions of 'v'.
"""
delm!(v::ENIgMaVert,idm) = setdiff!(v.make,idm);
"""
    delf!(v::ENIgMaVert,idf)

    Deletes 'idf' from the set of feed interactions of 'v'.
"""
delf!(v::ENIgMaVert,idf) = setdiff!(v.feed,idf);

function dele!(predator,predatorId,prey,preyId)   #deletes an eat interaction together with the reciprocal feed interaction
    dele!(predator,preyId);
    dele!(prey,predatorId);
end

function delm!(spec,specid,mod,modid)
    delm!(spec,modid);
    deln!(mod,specid);
end

Base.copy(vert::ENIgMaVert) = ENIgMaVert(copy(vert.eat),copy(vert.need),copy(vert.make),copy(vert.feed));

abstract type AbstractENIgMaGraph <: Graphs.AbstractGraph{Int} end

"""
    Defines a network in the ENIgMa-model.

# Extended help
    'vert::Dict{Int,ENIgMaVert}': Stores id=>vertex pairs of vertices in network.
    'spec::Set{Int}': Stores ids of species in the network.
    'basalRes::Set{Int}': Stores ids of basal resources in the network.
    'mods::Set{Int}': Stores ids of modifiers in the network.

    Provides 'hasv[id]' to efficiently check if vetex with 'id' is in the network.
    Provides 'hasspec[id]', 'hasBasalRes[id]' and 'hasMod[id]' to efficiently check if spec, basal resource and modifier with 'id' are in the network respectively.
    This can also be used to check if vertex with id is a species, a basal resource or a modifier.
"""
mutable struct ENIgMaGraph <: AbstractENIgMaGraph
    estsize::Int
    idmanager::IdManager
    hasv::BitVector             #is vert with id in Graph
    hasspec::BitVector          #is species with id in graph (is vert with id a species (in the graph))
    hasBasalRes::BitVector
    hasMod::BitVector
    spec::Set{Int}
    basalRes::Set{Int}
    mods::Set{Int}
    vert::Dict{Int,ENIgMaVert}  #stores vertices
    ENIgMaGraph(hasv::BitVector,hasspec::BitVector,hasBasalRes,hasMod,spec::Set{Int},basalRes,mods,vert::Dict{Int,ENIgMaVert}, idmanager::IdManager) = 
        new(length(hasv),idmanager,hasv,hasspec,hasBasalRes,hasMod,spec,basalRes,mods,vert);
end;

"""
    ENIgMaGraph(estsize::Int,idmanager::IdManager)

    Creates an empty Network wit estimated size 'estsize'.
"""
function ENIgMaGraph(estsize::Int,idmanager::IdManager)
    hasv = falses(estsize);
    hasspec = falses(estsize);
    hasBasalRes = falses(estsize)
    hasMod = falses(estsize)
    spec = Set{Int}();
    sizehint!(spec,estsize)
    basalRes = Set{Int}();
    sizehint!(basalRes,estsize)
    mods = Set{Int}();
    sizehint!(mods,estsize)
    vert = Dict{Int,ENIgMaVert}();
    sizehint!(vert,estsize);
    return ENIgMaGraph(hasv, hasspec, hasBasalRes, hasMod, spec, basalRes, mods, vert, idmanager);
end

ENIgMaGraph(estsize::Int) = ENIgMaGraph(estsize, IdManager(0));

getNextId!(g::ENIgMaGraph) = getNextId!(g.idmanager);

Base.iterate(g::ENIgMaGraph,i...) = iterate(g.vert,i...);
Base.length(g::ENIgMaGraph) = length(g.vert);
#Base.eltype(g::ENIgMaGraph) = eltype(g.vert);

#this allows to acces the vert dictionary with [] directly used on the network (poolnet[4] instead of poolnet.vert[4] etc)
Base.getindex(g::ENIgMaGraph,key...) = getindex(g.vert,key...);


#allows to compare two networks
function Base.:(==)(g1::ENIgMaGraph,g2::ENIgMaGraph)
    g1.idmanager.maxid != g2.idmanager.maxid && return false
    estSize = min(g1.estsize,g2.estsize)
    if g1.hasv[1:estSize] != g2.hasv[1:estSize] ||
       g1.hasspec[1:estSize] != g2.hasspec[1:estSize] ||
       g1.hasBasalRes[1:estSize] != g2.hasBasalRes[1:estSize] ||
       g1.hasMod[1:estSize] != g2.hasMod[1:estSize] ||
       g1.spec != g2.spec ||
       g1.mods != g2.mods ||
       g1.basalRes != g2.basalRes
        return false
    end
    for (id,v) in g1
        cv = g2[id]
        if v.eat != cv.eat || 
           v.need != cv.need ||
           v.make != cv.make ||
           v.feed != cv.feed2
            return false
        end
    end
    return true;
end

#called if the network needs to be enlarged as the network grew bigger than expected
function enlarge!(g::ENIgMaGraph,estsize::Int)
    estsize <= g.estsize && error("tried to shrink graph with enlarge!");    #failsave, shouldnt be neccessary
    g.hasv = vcat(g.hasv,falses(estsize - g.estsize));
    g.hasspec = vcat(g.hasspec,falses(estsize - g.estsize));
    g.hasBasalRes = vcat(g.hasBasalRes,falses(estsize - g.estsize));
    g.hasMod = vcat(g.hasMod,falses(estsize - g.estsize));
    g.estsize = estsize;
end

#returns number of species in network
getNumSpec(g::ENIgMaGraph) = length(g.spec);
getNumBasalRes(g::ENIgMaGraph) = length(g.basalRes)
getNumMods(g::ENIgMaGraph) = length(g.mods)

"""
    addspec!(g::ENIgMaGraph,id,v::ENIgMaVert)

    Adds species 'v' with id 'id' to 'g' as it is.

    Does not take care of further implications like feed interactions or creating objects. See colonize! for that.
"""
function addSpec!(g::ENIgMaGraph,id,v::ENIgMaVert)
    id > g.estsize && enlarge!(g,Int(round(id*enlargementfactor)))
    g.vert[id] = v;
    g.hasv[id] = true;
    push!(g.spec,id);
    g.hasspec[id] = true;
end

replacespec!(g::ENIgMaGraph,id,v::ENIgMaVert) = addSpec!(g::ENIgMaGraph,id,v::ENIgMaVert);     #here in case at some point erasing of old vert needs work might be optimized by having just g.vert[id] = v; would not make sure that vert has been there though

"""
    addBasalRes!(g::ENIgMaGraph,id,v::ENIgMaVert)

    Adds basal resource 'v' with id 'id' to 'g'.
"""
function addBasalRes!(g::ENIgMaGraph,id,v::ENIgMaVert)
    id > g.estsize && enlarge!(g,Int(round(id*enlargementfactor)))
    g.vert[id] = v;
    g.hasv[id] = true;
    push!(g.basalRes,id)
    g.hasBasalRes[id] = true;
end

"""
    addMod!(g::ENIgMaGraph,id,v::ENIgMaVert)

    Adds modifier/object 'v' with id 'id' to 'g'.
"""
function addMod!(g::ENIgMaGraph,id,v::ENIgMaVert)
    id > g.estsize && enlarge!(g,Int(round(id*enlargementfactor)))
    g.vert[id] = v;
    g.hasv[id] = true;
    push!(g.mods,id)
    g.hasMod[id] = true;
end

"""
    delv!(g::ENIgMaGraph,id::Int)

    Removes vertex v (species or modifier) from g.
    Deleting all references to it in the network (except needs).
"""
function delv!(g::ENIgMaGraph,id::Int)
    del_v = g[id]
    if g.hasspec[id]    #is spec?
        for eid in del_v.eat    #notify prey of absence
            delf!(g[eid],id)
        end
        for mid in del_v.make   # notify engineered mods of absence
            deln!(g[mid],id)
        end
        setdiff!(g.spec,id)
        g.hasspec[id] = false;  #could possibly be optimized by differentiating between spec and mod
    elseif g.hasMod[id]
        setdiff!(g.mods,id)
        g.hasMod[id] = false
    elseif g.hasMod
        setdiff!(g.basalRes,id)
        g.hasBasalRes[id] = false
        @warn "Deleted basal resource with id $(id)."
    end
    for fid in del_v.feed   #notify predators of absence
        dele!(g[fid],id)
    end
    delete!(g.vert,id);
    g.hasv[id] = false;
end

Base.in(id::Int,g::ENIgMaGraph) = g.hasv[id];
#hasv(g::ENIgMaGraph,id) = g.hasv[id]; hasnt been used might have been smarter...
include("GraphInterface.jl")
include("convertENIgMaGraphs.jl")
include("ENIgMaEvents.jl")

#calculates the competitive strength of all vertices
function calcstrength(poolnet::ENIgMaGraph, colnet::ENIgMaGraph,ce,cn,cm,cf)
    for specId in colnet.spec
        spec = colnet[specId]                  #might be optimizable by looping over vals directly
        spec.strength = -ce*length(poolnet[specId].eat) - cf*length(spec.feed) + cn*sum(colnet.hasv[collect(spec.need)]) + cm*length(spec.make); #could be optimized by using view and or finding better way to index by set
    end
end

"""
    getprimext(poolnet::ENIgMaGraph, colnet::ENIgMaGraph,ce,cn,cm,cf)

    Return a vector of ids of species that can go extinct by primary extinction in 'g'.

    # Arguments
    - 'ce::Float64': strength bonus per eat interaction.
    - 'cn::Float64': strength bonus per need interaction.
    - 'cn::Float64': strength bonus per make interaction.
    - 'cf::Float64': strength bonus per feed interaction.
"""
function getprimext(poolnet::ENIgMaGraph, colnet::ENIgMaGraph,ce,cn,cm,cf,secextspec)
    calcstrength(poolnet,colnet,ce,cn,cm,cf);
    
    #at first every spec could go extinct
    ext = copy(colnet.hasspec)  #safes if species goes extinct
    ext[secextspec] .= false;    #secondary extinction has higher priority    
    #allocate for later
    maxids = Int64[];
    maxstrength = typemin(Int);

    #now let every vert check which species are their strongest consumer(s) and delete those from the list
    for (fid, vert) in colnet
        isempty(vert.feed) && continue
        
        #find all strongest consumers
        maxstrength = typemin(Int);
        for eid in vert.feed
            strength = colnet[eid].strength;
            if strength > maxstrength           # could possibly be optimized by better algorithm (not self made multi argmax)
                maxstrength = strength;
                empty!(maxids);
                push!( maxids, eid );
            elseif maxstrength <= strength
                push!( maxids, eid);
            end
        end

        ext[maxids] .= false;   #delete strongest consumers from extinction list
    end

    return [id for id in eachindex(ext) if ext[id]]; #return vec of all species ids that could go extinct
end

"""
    getsecext!(g::ENIgMaGraph,extspec = Int[], extobj = Int[])

    Returns the two vectors of all species (via 'extspec') and modifier (via 'extobj') ids that could secondarily go extinct.
"""
function getsecext!(g::ENIgMaGraph,extspec = Int[], extobj = Int[])
    empty!(extspec);    #might be optimized by using other data structure or giving a sizehint
    empty!(extobj);
    for (id,vert) in g
        if g.hasspec[id]    #first treat species
            if isempty(vert.eat) || !all(g.hasv[collect(vert.need)]) # at least one eat and all need relationships must be established and 
                push!(extspec,id);
            end
        elseif g.hasMod[id]
            !any(g.hasv[collect(vert.need)]) && push!(extobj,id);   #might be optimized by finding better way to index by set.. same in if above
        end
    end
    return (extspec,extobj); #might be slightly optimized by returning nothing
end

"""
    getpotcolonizers!(poolnet::ENIgMaGraph,colnet::ENIgMaGraph,colonizers = Int[])
    
    Return a list (via 'colonizers') of species ids from the pool network that could colonize the colony network.
"""
function getpotcolonizers!(poolnet::ENIgMaGraph,colnet::ENIgMaGraph,colonizers = Int[])
    potcolonizers = setdiff(poolnet.spec,colnet.spec)
    empty!(colonizers);     #might be optimized by creating a vector from the Bitvec or giving a sizehint
    for id in potcolonizers #would be faster to first find all true values? optimize
        vert = poolnet[id];
        if all(colnet.hasv[collect(vert.need)]) && any(colnet.hasv[collect(vert.eat)])
            push!(colonizers,id)
        end
    end
    return colonizers;
end

#Introduces modifier/object from pool network to colony
function introduceMod!(poolnet::ENIgMaGraph,colnet::ENIgMaGraph,id_mod)
    feed = Set{Int}()
    for idf in poolnet[id_mod].feed    #...and add feeds manualy
        if colnet.hasv[idf]
            adde!(colnet[idf],id_mod);
            push!(feed,idf);
        end
    end
    addMod!(colnet, id_mod, ENIgMaVert(Set{Int}(),Set{Int}(),Set{Int}(),feed)) #add modifier to colony
end

"""
    function colonize!(poolnet::ENIgMaGraph,colnet::ENIgMaGraph,colonizerid,colonizer = copy(poolnet.vert[colonizerid]))
    
    Introduces species 'colonizer' with id 'colonizerid' from 'poolnet' to 'colnet'.
    
    Notifies other vertices of presence and introduces modifiers if necessary.
"""
function colonize!(poolnet::ENIgMaGraph,colnet::ENIgMaGraph,colonizerid,colonizer = nothing)
    if colonizer === nothing
        colonizer = copy(poolnet.vert[colonizerid])
    end
    #adjust eat and feed list
    for ide in colonizer.eat            #might be optimized by putting that into extra function
        if colnet.hasv[ide]
            addf!(colnet[ide],colonizerid);
        else
            dele!(colonizer, ide);       #might be optimized by adding existing links instead of deleting non existent ones
        end
    end

    for idf in colonizer.feed
        if colnet.hasv[idf]
            adde!(colnet[idf],colonizerid)
        else
            delf!(colonizer, idf);       #might be optimized by adding existing links instead of deleting non existent ones
        end
    end

    #add modifiers
    for id_mod in colonizer.make
        if !colnet.hasv[id_mod]  #if modifier doesnt exist yet create it
            introduceMod!(poolnet,colnet,id_mod)
        end
        addn!(colnet[id_mod],colonizerid)
    end

    addSpec!(colnet, colonizerid, colonizer);    # add species to colony
end



"""
    getinteractiontype(g:ENIgMaGraph, id1, id2)

    Return type of interaction (0 = ignore, 1 = eat, 2 = need, 3 = make) from 'id1' to 'id2' in 'g'.

    (feed interactions arent considered, as they go in the wrong direction)
"""
function getinteractiontype(g::ENIgMaGraph, id1, id2)
    v1 = g[id1]
    if id2 in v1.eat
        return 1;
    elseif id2 in v1.need
        return 2;
    elseif id2 in v1.make
        return 3;
    else
        return 0;
    end
end

#computes the difference in strength after a specific mutation
function getdeltastrength(old_int,new_int,change_in_interaction,ce,cn,cm,cpred)
    interaction_bonus = change_in_interaction ? Dict(0=>0,1=>-ce,2=>cn,3=>cm) : Dict(0=>0,1=>-cpred,2=>0,3=>0)
    return interaction_bonus[new_int] - interaction_bonus[old_int];
end


"""
    mutate!(poolnet::ENIgMaGraph, colnet::ENIgMaGraph, spmutid, intmutid, diverse, spints, obints, tallytable, evolutiontable, ce,cn,cm,cf)

    Realizes a random mutation of the species with id 'spmutid' in 'colnet' by changing one of its interactions with the vertex with id 'intmutid'.

    #Arguments
    - 'diverse': If true the mutated species is a completely new species. If false, it replaces the original if it is the stronger competitor and is discarded otherwise.
    - 'ce','cn','cm',cf': Factors that determine the competitive strength added by eat, need, make and feed interactions.

"""
function mutate!(poolnet::ENIgMaGraph, colnet::ENIgMaGraph, spmutid, intmutid, change_in_int, old_int, new_int, diverse, ce, cn, cm, cpred)
    if diverse == 1
        newid = getNextId!(poolnet);
        #if newid == 341
        #    println(c);
        #end
    else
        newid = spmutid
        if getdeltastrength(old_int,new_int,change_in_int,ce,cn,cm,cpred) < 0
            return newid;  #what would be the right tally for a discarded mutation?
        end
    end
    
    mutspecpool = copy(poolnet[spmutid]);
    interactorpool = poolnet[intmutid];
    mutspeccol = copy(colnet[spmutid]);
    interactorcol = colnet[intmutid];

    if change_in_int       
        #delete old interaction
        if old_int == 1
            dele!(mutspecpool, intmutid);
            dele!(mutspeccol, intmutid);
            if diverse == 0
                delf!(interactorpool, spmutid);
                delf!(interactorcol, spmutid);
            end
        elseif old_int == 2
            deln!(mutspecpool, intmutid);
            deln!(mutspeccol, intmutid);
        elseif old_int == 3
            delm!(mutspecpool, intmutid);
            delm!(mutspeccol, intmutid);
            if diverse == 0
                deln!(interactorpool, spmutid);
                deln!(interactorcol, spmutid);
            end
        end

        #add new interaction
        if new_int == 1
            adde!(mutspecpool, intmutid);
            adde!(mutspeccol, intmutid);
            addf!(interactorpool, newid);
            addf!(interactorcol, newid);
        elseif new_int == 2
            addn!(mutspecpool,intmutid);
            addn!(mutspeccol,intmutid);
        elseif new_int == 3
            addm!(mutspecpool, intmutid);
            addm!(mutspeccol, intmutid);
            addn!(interactorpool, newid);
            addn!(interactorcol, newid);
        end

    else #change out interaction
        #delete old interaction
        if old_int == 1
            if diverse == 0
                dele!(interactorpool, spmutid);
                dele!(interactorcol, spmutid);
            end
            delf!(mutspecpool, intmutid);
            delf!(mutspeccol, intmutid);
        elseif old_int == 2 && diverse == 0
            deln!(interactorpool, spmutid);
            deln!(interactorcol, spmutid);
        end

        #add new interaction
        if new_int == 1
            adde!(interactorpool, newid);
            adde!(interactorcol, newid);
            addf!(mutspecpool, intmutid);
            addf!(mutspeccol, intmutid);
        elseif new_int == 2
            addn!(interactorpool, newid);
            addn!(interactorcol, newid);
        end
    end

    if diverse == 0    #if diversification is disabled keep the mutated spec if its stronger then the original
        replacespec!(poolnet,spmutid,mutspecpool);
        replacespec!(colnet,spmutid,mutspeccol);
    elseif diverse == 1
        for eid in mutspecpool.feed         #add new id to other interactees of mutant spec
            adde!(poolnet[eid],newid)
        end
        for eid in mutspeccol.feed
            adde!(colnet[eid],newid)
        end
        for fid in mutspecpool.eat
            addf!(poolnet[fid],newid)
        end
        for fid in mutspeccol.eat
            addf!(colnet[fid],newid)
        end
        for mid in mutspecpool.make
            addn!(poolnet[mid],newid)
            addn!(colnet[mid],newid)
        end
        addSpec!(poolnet,newid,mutspecpool);
        addSpec!(colnet,newid,mutspeccol);
    end
    return newid;
end

include("setuppool.jl")
include("assemblyevo3.jl")

"""
    recreatecolnetdiverse(poolnet::ENIgMaGraph,it,ids,maxid,globextspec)

    Recreates the state of the colony network at itteration 'it' in the assembly process
        using the final pool network 'poolnet'.

    #Arguments:
    - 'ids::Bitvector': Stores which species are in the colony at the given itteration. 
        Is true for all ids of species in the colony network.
    - 'maxid::Int': The maximal id given to a vertex at the given itteration.
    - 'globextspec::Dict{Int,Pair{Int,ENIgMaVert}}': Stores all (id => vertex) pairs
         of globally extinct species using the itteration they went extinct as keys.
    """
function recreatecolnetdiverse(poolnet::ENIgMaGraph,it,ids,maxid::Int,globextspec)
    poolnet = deepcopy(poolnet)
    #make sure the right species are in the poolnet.
    for (extit,(id,v)) in globextspec   # Loop over all species extinct after simulation
        if extit <= it                  # already extinct at itteration that is recreated?
            if poolnet.hasv[id]         # is currently in poolnet? (eg due to earlier use of this function)
                delv!(poolnet,id);
            end
        elseif !poolnet.hasv[id]        # not yet extinct, but not in net?
            colonize!(poolnet,poolnet,id,v);
        end
    end
    colnet::ENIgMaGraph = ENIgMaGraph(poolnet.estsize,IdManager(maxid));    #create empty net
    for id in eachindex(ids)                #add vertices one by one
        if id == 265
            print("debug")
        end
        if ids[id]                          #is vert present?                                   
            if poolnet.hasspec[id]          #is it a spec?
                spec = copy(poolnet[id]);
                for nid in spec.need        #remove needs to verts, that don't exist at that point
                    if nid > maxid
                        deln!(spec, nid);
                    end
                end
                #use colonize! to introduce spec (it takes care of eats, feeds and makes and the respective modifier needs automatically)
                colonize!(poolnet,colnet,id,spec)   
            elseif poolnet.hasMod[id]                           #its a modifier
                if !colnet.hasv[id]         #only add modifiers, that haven't yet been introduced by a species
                    introducemod!(poolnet,colnet,id);   #introducemod! takes care of feeds
                end
            else
                basalRes = copy(poolnet[id])
                for fId in basalRes.feed
                    if !colnet.hasv[fId]
                        delf!(basalRes,fId)
                    end
                end
                addBasalRes!(colnet,id,basalRes)
            end
        end
    end
    return colnet;
end

recreatecolnetdiverse(poolnet,it,vertsInColony,maxids,globExtSpec) =
    recreatecolnetdiverse(poolnet,it,vertsInColony[:,it],maxids[it],globExtSpec)

recreatecolnetdiverse(simulationData::ENIgMaSimulationData,it) = 
    recreatecolnetdiverse(simulationData.poolnet,it,simulationData.vertsInColony,
        simulationData.maxids,simulationData.globExtSpec)

function getTrophicLevels_shortestPath(net::ENIgMaGraph)
    spec = copy(net.spec)
    trophic_lvls = []
    prev_lvl = net.basalRes
    curr_lvl = Set{Int}()
    while true
        for id in prev_lvl
            union!(curr_lvl,net[id].feed)
        end
        for lvl in trophic_lvls
            setdiff!(curr_lvl,lvl)
        end
        setdiff!(spec,curr_lvl)
        if isempty(curr_lvl)
            break
        end
        push!(trophic_lvls,copy(curr_lvl))
        curr_lvl, prev_lvl = empty!(prev_lvl),curr_lvl  #swap variables
    end
    return trophic_lvls
end

function getTrophicLevels(net)
    pathLengths = Dict{Int,Vector{Int}}()
    function followPath(path)
        currId = path[end]
        if haskey(pathLengths,currId)
            push!(pathLengths[currId], length(path))
        else
            pathLengths[currId] = [length(path)]
        end

        for specId in net[currId].feed
            if !(specId in path)
                followPath(push!(path, specId))
            end
        end
    end

    for basalResId in net.basalRes
        for specId in net[basalResId].feed
            followPath([specId])
        end
    end
    return Dict{Int,Float64}([specId => mean(lengths) for (specId,lengths) in pathLengths])
end

function getConnectedSpec(net)
    isConnected = falses(net.idmanager.maxid)
    function getConnected(id,isConnected,net)
        for fId in net[id].feed
            if !isConnected[fId]
                isConnected[fId] = true
                getConnected(fId,isConnected,net)
            end
        end
    end
    for basalResId in net.basalRes
        getConnected(basalResId,isConnected,net)
    end
    return [ind for ind in eachindex(isConnected) if isConnected[ind]]
end

function _getDegrees(g, linkTypes)
    degrees = zeros(Int,getNumSpec(g))
    for (ind,specId) in enumerate(g.spec)
        spec = g[specId]
        for type in linkTypes
            degrees[ind] += length(getfield(spec,type))
        end
    end
    return degrees
end

function getDegrees(g::ENIgMaGraph, linkType::Symbol)
    if linkType === :all || linkType === :in
        return _getDegrees(g, (:eat,:need,:make))
    end
    linkType == :pred && (linkType = :feed) 
    if linkType ∉ (:eat,:need,:make,:feed);
        error(":$(linkType) could not be interpreted as a type of interactions.\n\t Valid are: :all,:in,:eat,:need,:make,:feed and :pred.")
    end
    return _getDegrees(g,(linkType,))
end

function getDegreeDist(degrees::Vector{Int})
    uniqueDegrees = unique(degrees)
    maxDegree = maximum(uniqueDegrees)

    degreeCounts = zeros(Int,maxDegree+1)
    for degree in uniqueDegrees
        degreeCounts[degree + 1] = count(==(degree),degrees)
    end
    return 0:maxDegree,degreeCounts
end

getDegreeDist(g::ENIgMaGraph, linkType::Symbol) = getDegreeDist(getDegrees(g,linkType));

function _modifiedSDSimilarities(g,specIds,types)
    S = length(specIds)
    SDs = zeros(Int(S*(S-1)/2))
    for (it1,id1) in enumerate(specIds[1:end-1])
        for (it2,id2) in enumerate(specIds[it1+1:end])
            spec1 = g[id1]
            spec2 = g[id2]
            #S*(S-1)/2 - (S - it1 + 1)(S - it1)/2 + it2 ==
            #(S^2 - S - S^2 + S*it1 + S*it1 - it1^2 - S + it1)/2 + it2  ==
            #S*it1 - S - it1*(it1 - 1)/2 + it2 (index: big triangle - small triangle (see as linear indices in lower triangular matrix))
            SDs[S*it1 - S - it1*(it1 - 1)÷2 + it2] = sum(2*(length(intersect(getfield(spec1,type),getfield(spec2,type)))) for type in types)/
                    sum(length(getfield(spec1,type)) + length(getfield(spec2,type)) for type in types)
        end
    end
    return SDs,specIds
end

"""
modifiedSDSimilarities(g::ENIgMaGraph, specIds, linkType::Symbol)

Computes a modified version of the Soerensen-Dice Similarity Coefficient for all distinct
pairs of distinct species.

This modified Soerensen-Dice Coefficient of two species is defined to be two times the 
number of links with type defined by 'linkType' they have in common
(meaning the same type of link to the same node) divided by the sum of their total
number of links with types defined by 'linkTypes'.

"""
function modifiedSDSimilarities(g::ENIgMaGraph, specIds, linkType::Symbol)
    if linkType === :all || linkType === :in
        return _modifiedSDSimilarities(g, specIds, (:eat,:need,:make))
    end
    linkType == :pred && (linkType = :feed) 
    if linkType ∉ (:eat,:need,:make,:feed);
        error(":$(linkType) could not be interpreted as a type of interactions.\n\t Valid are: :all,:in,:eat,:need,:make,:feed and :pred.")
    end
    return _modifiedSDSimilarities(g,specIds,(linkType,))
end

modifiedSDSimilarities(g::ENIgMaGraph) = 
    modifiedSDSimilarities(g,collect(g.spec), :all)

modifiedSDSimilarities(g::ENIgMaGraph,specIds) =
    modifiedSDSimilarities(g::ENIgMaGraph,specIds, :all)

modifiedSDSimilarities(g::ENIgMaGraph, linkType::Symbol) =
     modifiedSDSimilarities(g,collect(g.spec),linkType)



end #module