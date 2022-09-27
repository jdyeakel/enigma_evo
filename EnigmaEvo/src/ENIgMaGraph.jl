# defines a Graph type with multiple edge types stored in seperate adjacency lists (Sets actually)

#optimization ideas: include strength calculations in add and del calls

#performance questions: is there a good graph representation out there? Sets good choice in vert? dict good choice in graph

module ENIgMaGraphs

export ENIgMaGraph, ENIgMaVert, IdManager       #might be put in begin block? kinda didnt work
export adde!, addn!, addf!, addm!, dele!, deln!, delm!, delf!
export numspec
export addspec!,replacespec!,addmod!, delv!
export getprimext, getsecext!, getpotcolonizers!
export colonize!,mutate!
export getnextid!

const enlargementfactor = 1.1; #controls how much buffer is added if estsize has to be increased

mutable struct IdManager
    maxid::Int

    IdManager(offset) = new(offset)
end;

setmaxid!(idmanager,maxid) = idmanager.maxid = maxid;

function getnextid!(idmanager::IdManager) #optimization idea: have IdManager save ids of globally extinct species to reassign them?
    idmanager.maxid += 1;
    return idmanager.maxid;
end
mutable struct ENIgMaVert
    eat::Set{Int}
    need::Set{Int}
    make::Set{Int}
    feed::Set{Int}

    strength::Float64

    ENIgMaVert() = new(Set{Int}(),Set{Int}(),Set{Int}(),Set{Int}(),NaN)
    ENIgMaVert(eat,need,make,feed) = new(eat,need,make,feed,NaN)
end;

#if changed: these have not been used all the time
adde!(v::ENIgMaVert,ide) = push!(v.eat,ide);
addn!(v::ENIgMaVert,idn) = push!(v.need,idn);
addm!(v::ENIgMaVert,idm) = push!(v.make,idm);
function addf!(v::ENIgMaVert,idf)
    #if idf == 341
    #    println("debug addf");
    #end
     push!(v.feed,idf); 
end

dele!(v::ENIgMaVert,ide) = setdiff!(v.eat,ide);
deln!(v::ENIgMaVert,idn) = setdiff!(v.need,idn);
delm!(v::ENIgMaVert,idm) = setdiff!(v.make,idm);
delf!(v::ENIgMaVert,idf) = setdiff!(v.feed,idf);

Base.copy(vert::ENIgMaVert) = ENIgMaVert(copy(vert.eat),copy(vert.need),copy(vert.make),copy(vert.feed));

mutable struct ENIgMaGraph
    estsize::Int
    idmanager::IdManager
    hasv::BitVector             #is vert with id in Graph
    hasspec::BitVector          #is species with id in graph (is vert with id a species (in the graph))
    spec::Set{Int}
    vert::Dict{Int,ENIgMaVert}  #stores verts
    ENIgMaGraph(hasv::BitVector,hasspec::BitVector,spec::Set{Int},vert::Dict{Int,ENIgMaVert}, idmanager::IdManager) = new(length(hasv),idmanager,hasv,hasspec,spec,vert);
end;

function ENIgMaGraph(estsize::Int,idmanager::IdManager)
    hasv = falses(estsize);
    hasspec = falses(estsize);
    spec = Set{Int}();
    sizehint!(spec,estsize)
    vert = Dict{Int,ENIgMaVert}();
    sizehint!(vert,estsize);
    return ENIgMaGraph(hasv,hasspec,spec,vert,idmanager);
end

ENIgMaGraph(estsize::Int) = ENIgMaGraph(estsize, IdManager(0));

getnextid!(g::ENIgMaGraph) = getnextid!(g.idmanager);

Base.iterate(g::ENIgMaGraph,i...) = iterate(g.vert,i...);
Base.length(g::ENIgMaGraph) = length(g.vert);
Base.eltype(g::ENIgMaGraph) = eltype(g.vert);

Base.getindex(g::ENIgMaGraph,key...) = getindex(g.vert,key...);

function enlarge!(g::ENIgMaGraph,estsize::Int)
    estsize <= g.estsize && error("tried to shrink graph with enlarge!");    #failsave, shouldnt be neccessary
    g.hasv = vcat(g.hasv,falses(estsize - g.estsize));
    g.hasspec = vcat(g.hasspec,falses(estsize - g.estsize));
    g.estsize = estsize;
end

#returns number of species in network
numspec(g::ENIgMaGraph) = length(g.spec);#sum(g.hasspec);

#add species to graph
function addspec!(g::ENIgMaGraph,id,v::ENIgMaVert)
    id > g.estsize && enlarge!(g,Int(round(id*enlargementfactor)))     #failsave, shouldnt be neccessary
    g.vert[id] = v;
    g.hasv[id] = true;
    push!(g.spec,id);
    g.hasspec[id] = true;
end

replacespec!(g::ENIgMaGraph,id,v::ENIgMaVert) = addspec!(g::ENIgMaGraph,id,v::ENIgMaVert);     #here in case at some point erasing of old vert needs work might be optimized by having just g.vert[id] = v; would not make sure that vert has been there though

#add abiotic modifier to graph
function addmod!(g::ENIgMaGraph,id,v::ENIgMaVert)
    id > g.estsize && enlarge!(g,Int(round(id*enlargementfactor)))     #failsave, shouldnt be neccessary
    g.vert[id] = v;
    g.hasv[id] = true;
end

function delv!(g::ENIgMaGraph,id::Int)
    del_v = g[id]
    for fid in del_v.feed
        dele!(g[fid],id)
    end
    for eid in del_v.eat
        delf!(g[eid],id)
    end
    delete!(g.vert,id);
    g.hasv[id] = false;
    if g.hasspec[id]
        setdiff!(g.spec,id)
        g.hasspec[id] = false;  #could possibly be optimized by differentiating between spec and mod
    end
end

Base.in(id::Int,g::ENIgMaGraph) = g.hasv[id];
#hasv(g::ENIgMaGraph,id) = g.hasv[id]; hasnt been used might have been smarter...

#calculates the competitive strength of all vertices
function calcstrength(g::ENIgMaGraph,ce,cn,cm,cf)
    for (_,vert::ENIgMaVert) in g                   #might be optimizable by looping over vals directly
        vert.strength = -ce*length(vert.eat) - cf*length(vert.feed) + cn*sum(g.hasv[collect(vert.need)]) + cm*length(vert.make); #could be optimized by using view and or finding better way to index by set
    end
end

"""
    getprimext(g::ENIgMaGraph,ce,cn,cm,cf)

    Return a vector of ids of species that can go extinct by primary extinction in 'g'.

    # Arguments
    - 'ce::Float64': strength bonus per eat interaction.
    - 'cn::Float64': strength bonus per need interaction.
    - 'cn::Float64': strength bonus per make interaction.
    - 'cf::Float64': strength bonus per feed interaction.
"""
function getprimext(g::ENIgMaGraph,ce,cn,cm,cf)
    calcstrength(g,ce,cn,cm,cf)
    
    #at first every vert could go extinct
    ext = copy(g.hasspec)   #safes if species goes extinct

    #allocate for later
    maxids = Int64[];
    maxstrength = typemin(Int);

    #now let every vert check which species are their strongest consumer(s) and delete those from the list
    for (fid, vert) in g
        isempty(vert.feed) && continue

        # there is no competition over basal resource by model definition
        if fid == 1        
            for eid in vert.feed
                ext[eid] = false;
            end
            continue
        end
        
        #find all strongest consumers
        maxstrength = typemin(Int);
        for eid in vert.feed
            #if !haskey(g.vert,eid)#for debugging
            #    erorr()
            #end
            strength = g[eid].strength;
            if strength > maxstrength           # could possibly be optimized by better algorithm (not self made multi argmax)
                maxstrength = strength;
                empty!(maxids);
                push!( maxids, eid );
            elseif maxstrength >= strength
                push!( maxids, eid);
            end
        end

        ext[maxids] .= false;   #delete strongest consumers from extinction list
    end

    return [id for id = 1:length(ext) if ext[id]]; #return vec of all species ids that could go extinct
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
        else
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
    potcolonizers = xor.(poolnet.hasspec, colnet.hasspec);   #xor is not exactly what we want. If we do everything right that shouldnt matter though
    empty!(colonizers);     #might be optimized by creating a vector from the Bitvec or giving a sizehint
    for id in eachindex(potcolonizers) #would be faster to first find all true values? optimize
        if potcolonizers[id]
            vert = poolnet[id];
            if all(colnet.hasv[collect(vert.need)]) && any(colnet.hasv[collect(vert.eat)])
                push!(colonizers,id)
            end
        end
    end
    return colonizers;
end

function colonize!(poolnet::ENIgMaGraph,colnet::ENIgMaGraph,colonizerid)
    colonizer = copy(poolnet[colonizerid])

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
        if colnet.hasv[id_mod]  #if modifier already there just add need to modifier
            addn!(colnet[id_mod],colonizerid)
        else    #otherwise create a new one...
            feed = Set{Int}()
            for idf in poolnet[id_mod].feed    #...and add feeds manualy
                if colnet.hasv[idf]
                    adde!(colnet[idf],id_mod);
                    push!(feed,idf);
                end
            end
            addmod!(colnet, id_mod, ENIgMaVert(Set{Int}(),Set{Int}(colonizerid),Set{Int}(),feed)) #add modifier to colony
        end
    end

    addspec!(colnet, colonizerid, colonizer);    # add species to colony
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


function mutate!(poolnet::ENIgMaGraph, colnet::ENIgMaGraph, spmutid, intmutid, diverse, spints, obints, tallytable, evolutiontable, ce,cn,cm,cf)
    mutspecpool = copy(poolnet[spmutid]);
    interactorpool = poolnet[intmutid];
    mutspeccol = copy(colnet[spmutid]);
    interactorcol = colnet[intmutid];

    if diverse == 1
        newid = getnextid!(poolnet);
        #if newid == 341
        #    println("debug");
        #end
    else
        newid = spmutid
    end
    
    #get current interaction type (faster in colnet)
    oldint_in = getinteractiontype(colnet, spmutid, intmutid);
    
    #how intmut interacts with spmut
    oldint_out = getinteractiontype(colnet, intmutid, spmutid);

    #intm_mut = copy(intm);

    deltastrength = 0;

    #is the interacting vertex a species?
    if colnet.hasspec[intmutid]
        newint_in = rand(setdiff(spints,oldint_in));     #could be optimized by choosing after in and out is chosen
        newint_out = rand(setdiff(spints,oldint_out));
        
        # Select evolution of either in degree interaction or out-degree interaction
        evol_type_draw = rand();

        #every case done by hand might be the fastest way though (if debugging is painfull it might be easier to just exchange indcs for different cases...guess just use symetry)
        if evol_type_draw < 0.5 #change in degree?
            newint = newint_in;
            oldint = oldint_in;
            
            #delete old interaction
            if oldint == 1
                dele!(mutspecpool, intmutid);
                dele!(mutspeccol, intmutid);
                deltastrength += ce;
                if diverse == 0
                    delf!(interactorpool, spmutid);
                    delf!(interactorcol, spmutid);
                end
            elseif oldint == 2
                deln!(mutspecpool, intmutid);
                deln!(mutspeccol, intmutid);
                deltastrength -= cn;
            end

            #add new interaction
            if newint == 1
                adde!(mutspecpool, intmutid);
                adde!(mutspeccol, intmutid);
                deltastrength -= ce;
                addf!(interactorpool, newid);
                addf!(interactorcol, newid);
            elseif newint == 2
                addn!(mutspecpool,intmutid);
                addn!(mutspeccol,intmutid);
                deltastrength += cn;
            end

            #intm_mut[spmut,intmut] = newint;
            evol_type = 1.;
        else #change out degree
            newint = newint_out;
            oldint = oldint_out;
            
            #delete old interaction
            if oldint == 1
                if diverse == 0
                    dele!(interactorpool, spmutid);
                    dele!(interactorcol, spmutid);
                end
                delf!(mutspecpool, intmutid);
                delf!(mutspeccol, intmutid);
                deltastrength += cf;
            elseif oldint == 2 && diverse == 0
                deln!(interactorpool, spmutid);
                deln!(interactorcol, spmutid);
            end

            #add new interaction
            if newint == 1
                adde!(interactorpool, newid);
                adde!(interactorcol, newid);
                addf!(mutspecpool, intmutid);
                addf!(mutspeccol, intmutid);
                deltastrength -= cf;
            elseif newint == 2
                addn!(interactorpool, newid);
                addn!(interactorcol, newid);
            end

            evol_type = 2.;
        end
    #If the mutated interaction is with an object
    else
        #For objects, randomly choose ignore, eat, need, make
        newint_in = rand(setdiff(obints,oldint_in))[1];
        newint = newint_in;
        oldint = oldint_in;
        
        #delete old interaction
        if oldint == 1
            dele!(mutspecpool, intmutid);
            dele!(mutspeccol, intmutid);
            deltastrength += ce;
            if diverse == 0
                delf!(interactorpool, spmutid);
                delf!(interactorcol, spmutid);
            end
        elseif oldint == 2
            deln!(mutspecpool, intmutid);
            deln!(mutspeccol, intmutid);
            deltastrength -= cn;
        elseif oldint == 3
            delm!(mutspecpool, intmutid);
            delm!(mutspeccol, intmutid);
            deltastrength -= cm;
            if diverse == 0
                deln!(interactorpool, spmutid);
                deln!(interactorcol, spmutid);
            end
        end
        
        #add new interaction
        if newint == 1
            adde!(mutspecpool, intmutid);
            adde!(mutspeccol, intmutid);
            deltastrength -= ce;
            addf!(interactorpool, newid);
            addf!(interactorcol, newid);
        elseif newint == 2
            addn!(mutspecpool,intmutid);
            addn!(mutspeccol,intmutid);
            deltastrength += cn;
        elseif newint == 3
            addm!(mutspecpool, intmutid);
            addm!(mutspeccol, intmutid);
            deltastrength += cm;
            addn!(interactorpool, newid);
            addn!(interactorcol, newid);
        end            
        evol_type = 1.;
    end

    # ACCEPT REGARDLESS
    #ebmut,nbmut,nb0mut,mbmut = intbool(intm_mut);
    tally = tallytable[(evolutiontable[1,:] .== oldint) .& (evolutiontable[2,:] .== newint)][1];

    #Record evolution type (in degree vs. out degree) on backend of tally
    tally += evol_type*0.01;

    if diverse == 0 && deltastrength > 0
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
        addspec!(poolnet,newid,mutspecpool);
        addspec!(colnet,newid,mutspeccol);
    end

    return tally;#(intm_mut, ebmut, nbmut, nb0mut, mbmut, tally)
end

end #module
#=
function test()
    g = ENIgMaGraph(10);
    for i = 1:8
        addspec!(g,i,ENIgMaVert());
    end
    println(g[2].eat)
    for (_,v) in g
        println(v);
    end
end
=#