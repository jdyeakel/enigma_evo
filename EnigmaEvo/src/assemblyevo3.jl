abstract type ENIgMaSimulationData end

"""
    Data type that stores the generated data of a simulation of the ENIgMaModel.

# Extended help
    #Fields
    -'poolnet::ENIgMaGraph': The pool network at the end of the simulation.
    -'colnet::ENIgMaGraph': The colony network at the end of the simulation.
    -'phyloTree::ManyRootTree' The phylogenetic tree generated if diversification is on.
    -'specRich::Vector{Int}': The colony's species richness in each itteration.
    -'rich::Vector{Int}': The number of vertices in the colony for each itteration.
    -'pool::Vector{Int}': The number of vertices in the pool for each itteration.
    -'mstrength::Vector{Float64}': The (somewhat normalized) mean strength for each itteration.
    -'clock::Vector{Float64}': The time since the start of the simulation for each itteration.
    -'vertsInColony::BitArray': Stores which vertices are in the colony for each itteration.
        vertsInColony[id,it] is true if vertex with id 'id' is in colony at itteration 'it'.
    -'maxids::Vector{Int}': Saves the maximal vertex id for each itteration.
    -'globExtSpec::Dict{Int,Pair{Int,ENIgMaVert}}': Stores all (id => vertex) pairs
        of globally extinct species using the itteration they went extinct as keys.
    -'meanEats::Vector{Float64}': The mean out-degree of the eat-subnetwork of the colony.
    -'meanNeeds::Vector{Float64}': The mean out-degree of the need-subnetwork of the colony.
    -'meanEats_pool::Vector{Float64}': The mean out-degree of the eat-subnetwork of the pool.
    -'meanNeeds_pool::Vector{Float64}': The mean out-degree of the need-subnetwork of the pool.
    -'meanSpecEats::Vector{Float64}': The mean out-degree of the species-species eat-subnetwork of the colony.
    -'meanSpecNeeds::Vector{Float64}': The mean out-degree of the species-species need-subnetwork of the colony. 
    -'meanSpecEats_pool::Vector{Float64}': The mean out-degree of the species-species eat-subnetwork of the pool.
    -'meanSpecNeeds_pool::Vector{Float64}': The mean out-degree of the species-species need-subnetwork of the pool. 
    -'events::Vector{AbstractENIgMaEvent}': Stores the event of each itteration.
    -'trophLevels::Vector{Vector{Float64}}': The trophic levels in the colony for each itteration.
    -'nPrimExtSpec::Vector{Float64}': The number of species that could go primarily extinct for each itteration.
    -'nSecExtSpec::Vector{Float64}': The number of species that could go secondary extinct for each itteration.
    -'nColonizers::Vector{Float64}': The number of species that could colonize for each itteration.
"""
struct ENIgMaSimulationData_v3 <: ENIgMaSimulationData
    poolnet::ENIgMaGraph
    colnet::ENIgMaGraph
    phyloTree::ManyRootTree
    specRich::Vector{Int}
    rich::Vector{Int}
    pool::Vector{Int}
    mstrength::Vector{Float64}
    clock::Vector{Float64}
    vertsInColony::BitArray
    maxids::Vector{Int}
    globExtSpec::Dict{Int,Pair{Int,ENIgMaVert}}
    meanEats::Vector{Float64}
    meanNeeds::Vector{Float64}
    meanEats_pool::Vector{Float64}
    meanNeeds_pool::Vector{Float64}
    meanSpecEats::Vector{Float64}
    meanSpecNeeds::Vector{Float64}
    meanSpecEats_pool::Vector{Float64}
    meanSpecNeeds_pool::Vector{Float64}
    events::Vector{AbstractENIgMaEvent}
    trophLevels::Vector{Vector{Float64}}
    nPrimExtSpec::Vector{Float64}
    nSecExtSpec::Vector{Float64}
    nColonizers::Vector{Float64}
end

function assemblyevo(poolnet::ENIgMaGraph, rates, maxits, cm, cn, ce, cf, diverse, restrict_colonization::Bool, initColNet = missing; createLog = true)
    #Is initial colony given?
    if initColNet === missing
        #create the colony network, use idmanager of pool network (no copy, a reference to the original)
        colnet::ENIgMaGraph = ENIgMaGraph(poolnet.estsize,poolnet.idmanager);
        #add basal resources as they are always there
        for basalResId in poolnet.basalRes
            addBasalRes!( colnet, basalResId, ENIgMaVert());
        end
    else
        colnet = initColNet
    end
    
    # set up some monitoring buffers
    specRich = zeros(Int64,maxits);
    rich = zeros(Int64,maxits);
    pool = zeros(Int64,maxits);
    clock = zeros(Float64,maxits);
    mstrength = zeros(Float64,maxits);
    events = Array{AbstractENIgMaEvent}(undef,maxits);
    if createLog
        vertsInColony = falses(poolnet.estsize,maxits);
        maxids = zeros(Int,maxits);
        globextspec = Dict{Int,Pair{Int,ENIgMaVert}}()    #Stores globally extinct species together with their ids. The keys are the iterations the species went extinct
        trophLevels = Vector{Vector{Float64}}()
    end

    #initialize phylogenetic tree
    phyloTree = ManyRootTree()  #create phylo genetic tree
    if diverse == 1.5  #in non diverse case phylogeny kinda weird
        ##create root nodes
        #createnodes!(phyloTree,Dict{String,Dict{String,Any}}("$(id)_root" => Dict{String,Any}("timestamp" => 0) for id in poolnet.spec))
        superRoot = createnode!(phyloTree,"superRoot", data=Dict{String,Any}("timestamp"=>0.0,"evolution" => ""))
        #create leaf nodes those will be kept and move along with time, see note in evlution passage
        createnodes!(phyloTree,Dict{String,Dict{String,Any}}("$id" => Dict{String,Any}("parentName"=>"superRoot", "version"=>1, "heritage"=>"") for id in poolnet.spec))
    end

    meanEats = zeros(Float64,maxits);
    meanNeeds = zeros(Float64,maxits);
    meanEats_pool = zeros(Float64,maxits);
    meanNeeds_pool = zeros(Float64,maxits);
    meanSpecEats = zeros(Float64,maxits);
    meanSpecNeeds = zeros(Float64,maxits);
    meanSpecEats_pool = zeros(Float64,maxits);
    meanSpecNeeds_pool = zeros(Float64,maxits);
    nPrimExtSpec = zeros(Float64,maxits)
    nSecExtSpec = zeros(Float64,maxits)
    nColonizers = zeros(Float64,maxits)

    #oldColonies = Dict{Int,ENIgMaGraph}()
    #oldPools = Dict{Int,ENIgMaGraph}()

    #mutstep = Float64[]#zeros(Float64,maxits);
    #evolvedstrength = Array{Float64}(undef,0);

    #evolutiontable = [[0 0 0 1 1 1 2 2 2 3 3 3];[1 2 3 0 2 3 0 1 3 0 1 2]];
    #tallytable = [4.1 4.2 5.1 4.3 4.4 5.2 4.5 4.6 5.3 6.1 6.2 6.3];
    interactionTable = Dict{Int, String}(0=>"i",1=>"e",2=>"n",3=>"m")

    # If diversification is turned off, rates.rext -> 0 do that by hand
    if diverse == 0
        rates = (rc = rates.rc,rsecext = rates.rsecext, rprimext = rates.rprimext, reo = rates.reo,revo = rates.revo,rext = 0.);
    end

    colonizers = Int[]; #later used to store ids of potential colonizers
    secextspec = Int[]; #later used to store ids of spec that could go extinct secondarily
    extobj = Int[];     #later used to store ids of objects that could go extinct

    t=0;
    it = 0;
    while it < maxits

        minstrength = ce*Float64(length(poolnet) - 1);    #assuming ce > cf, eat everything but yourself
        maxstrength = cm*Float64(getNumMods(colnet)) + cn*Float64(getNumSpec(colnet) + getNumBasalRes(colnet) - 1) - ce*1 - cf*0;  #assuming cm > cn

        #COUNT POTENTIAL COLONIZERS
        if restrict_colonization
            getpotcolonizers!(poolnet,colnet,colonizers);   #saves pot colonizers ids in colonizers
        else
            colonizers = setdiff(poolnet.spec,colnet.spec);
        end

        #COUNT SPECIES THAT COULD POTENTIALLY GO EXTINCT
        
        #1) SECONDARY EXTINCTIONS By not fulfilling Eat/Need thresholds
        getsecext!(colnet,secextspec,extobj);

        #2) PRIMARY EXTINCTIONS by competition
        primextspec = getprimext(poolnet,colnet,ce,cn,cm,cf,secextspec);

        #COUNT POTENTIAL evolutionary events (size of community)
        levo = max(getNumSpec(colnet) - 1,0); #means evolution can only occur if community has more than 1 species

        #COUNT POTENTIAL global extinction events (size of pool)
        lext =  getNumSpec(poolnet) - getNumSpec(colnet) #- (colnet.hasv[2] ? 0 : 1) #species 2 can't go exinct

        # Calculate the full Rate
        Rate = rates.rc*length(colonizers) + rates.rprimext*length(primextspec) + rates.rsecext*length(secextspec) + rates.reo*length(extobj) + rates.revo*levo + rates.rext*lext;
        if Rate == 0
            println("Warning: No action was taken at iteration $it as no action was possible.")
            println("         Maybe the system evolved to a state where no species can colonize and almost all went extinct?")
            println("         Aborting simulation prematurely.")
            break
        end

        # Calculate event probabilities
        probc = (rates.rc*length(colonizers))/Rate;
        probprimext = (rates.rprimext*length(primextspec))/Rate;
        probsecext = (rates.rsecext*length(secextspec))/Rate;
        probeo = (rates.reo*length(extobj))/Rate;
        probevo = (rates.revo*levo)/Rate;
        #probext = (rates.rext*lext)/Rate; not needed right now
       
        dt = 1/Rate;
        t += dt;
        it += 1;
        clock[it] = t;

        #Choose a random event
        re = rand();

        #DRAW COLONIZATION
        if re <= probc #(lcol/levents)

            #COLONIZATION FUNCTION
            idc = rand(colonizers);
            colonize!(poolnet,colnet,idc)

            events[it] = ColonizationEvent(idc)

        #DRAW primary EXTINCTION
        elseif re <= (probc + probprimext)
            #select species to go extinct
            sp_bye = rand(primextspec);
            delv!( colnet, sp_bye )

            #These are extinctions from competitive exclusion
            
            events[it] = PrimaryExtinctionEvent(sp_bye)

        #DRAW SECONDARY EXTINCTION
        elseif re <= (probc + probprimext + probsecext)

            #select species to go extinct
            sp_bye = rand(secextspec);
            delv!( colnet, sp_bye )

            #These are extinctions from structural failings
            events[it] = SecondaryExtinctionEvent(sp_bye)

        #DRAW OBJECT EXTINCTION
        elseif re <= (probc + probprimext + probsecext + probeo) #((lcol + lspext)/levents) && re < ((lcol + lspext + lobext)/levents)

            extObjId = rand(extobj)
            delv!( colnet, extObjId )

            events[it] = ObjectExtinctionEvent(extObjId)
        
        #DRAW EVOLUTION
        elseif re <= (probc + probprimext + probsecext + probeo + probevo)

            # EVOLUTION OF SPECIES IN THE COMMUNITY
            # SPECIES MUTATION
            # select FROM community
            spmutid = rand(colnet.spec);    #might be optimized, by finding smarter way, to draw from Set (a#same below for intmutid)
            #OR Mutate realized interactions
            intmutid = rand(collect(setdiff(keys(colnet.vert),spmutid)));    #could possibly be optimized by drawing in loop until valid value is drawn

            basalResInts = [0,1]
            spec_ints = [0,1,2];
            object_ints = [0,1,2,3];

            if colnet.hasspec[intmutid] #is interactee a spec?
                change_in_int = rand([true,false]); #change incoming interaction?
                if change_in_int
                    old_int = ENIgMaGraphs.getinteractiontype(colnet, spmutid, intmutid);
                else
                    old_int = ENIgMaGraphs.getinteractiontype(colnet, intmutid, spmutid);
                end
                new_int = rand(setdiff(spec_ints,old_int));
            else    #interactee is a modifier
                change_in_int = true;
                old_int = ENIgMaGraphs.getinteractiontype(colnet, spmutid, intmutid);
                if colnet.hasBasalRes[intmutid]
                    new_int = rand(setdiff(basalResInts,old_int)); #cant make primal resource (sun)
                else
                    new_int = rand(setdiff(object_ints,old_int));
                end
            end

            newId = mutate!(poolnet, colnet, spmutid, intmutid, change_in_int, old_int, new_int, diverse, ce, cn, cm, cf);
            


            if createLog
                if poolnet.estsize > size(vertsInColony)[1];      #has the network grown bigger than our buffer?
                    vertsInColony = vcat(vertsInColony,falses(poolnet.estsize - size(vertsInColony)[1],maxits));  #then extend the buffer
                end
            end

            events[it] = MutationEvent(spmutid,intmutid,newId,InteractionType(old_int),InteractionType(new_int), change_in_int);

            #update phylogenic Tree (approach a bit unintuitive.
            # I keep the leaf nodes unattached, just move them forward and let them save their current version
            # and where they are currently attached without actually attaching them)
            if diverse == 1.5
                leafNode = getnode(phyloTree,"$(spmutid)")      #get the leaf node of the mutating species to move it to the present        
                soonGrandParent = getnode(phyloTree,leafNode.data["parentName"])   #get its parent's and soon to be grandparent's name 

                #The plan: sgp----x  -->   sgp----np----x      (x = leafNode gets moved to present and stays leaf)
                #                                  `----y
                #create the node of the branching event, the new parent of the leafNode x
                newParentName = "$(spmutid)v$(leafNode.data["version"])"
                newParent = createnode!(phyloTree,newParentName,data=Dict{String,Any}("timestamp"=>t,"evolution"=>"$(change_in_int ? "i" : "o"):$(interactionTable[old_int])$(interactionTable[new_int])"))
                createbranch!(phyloTree,soonGrandParent,newParent,t - soonGrandParent.data["timestamp"])    #connect grandparent to parent (last argument length of branch) 

                #update leaf node and create new leaf for the evolved species
                createnode!(phyloTree,"$newId",data=Dict{String,Any}("parentName" => newParentName, "version" => 1, "heritage" => "$(leafNode.data["heritage"])_$(newParentName)"))
                leafNode.data["parentName"] = newParentName;
                leafNode.data["version"] += 1;
            end

        # Global EXTINCTION
        elseif re <= 1
            #Choose species subject to extinction
            #Do not choose species 2, which is protected
            globextid = rand(setdiff(poolnet.spec,colnet.spec));
            
            if colnet.hasv[globextid]
                delv!(colnet,globextid)
            end

            if createLog
                globextspec[it] = Pair(globextid,poolnet[globextid]);
            end

            delv!(poolnet,globextid);

            events[it] = GlobalExtinctionEvent(globextid) 

            #update Phylogenetic tree
            if diverse == 1.5
                nowExtNode = getnode(phyloTree,"$(globextid)")
                parent = getnode(phyloTree,nowExtNode.data["parentName"])
                createbranch!(phyloTree,parent,nowExtNode,t - parent.data["timestamp"])
                nowExtNode.data["extinct"] = true
                nowExtNode.data["timestamp"] = t
            end
        end
        
        #SAVE RESULTS IN BUFFERS
        numSpec = getNumSpec(colnet)
        numSpec_pool = getNumSpec(poolnet)

        if isempty(colnet.spec)
            meanEats[it] = meanNeeds[it] = meanSpecEats[it] = meanSpecNeeds[it] = mstrength[it] = NaN64
        else
            meanEats[it] = sum(length(colnet[id].eat) for id in colnet.spec)/numSpec;     #bit afraid of overflows...
            meanNeeds[it] = sum(length(colnet[id].need) for id in colnet.spec)/numSpec;
            meanSpecEats[it] = sum(count(colnet.hasspec[collect(colnet[id].eat)]) for id in colnet.spec)/numSpec;     #bit afraid of overflows...
            meanSpecNeeds[it] = sum(count(poolnet.hasspec[collect(colnet[id].need)]) for id in colnet.spec)/numSpec;
            mstrength[it] = (mean( v.strength for (id,v) in colnet if colnet.hasspec[id] ) - minstrength)/(maxstrength - minstrength); #make strengths positive with a minimum of 1
        end
        if isempty(poolnet.spec)
            meanEats_pool[it] = meanNeeds_pool[it] = meanSpecEats_pool[it] = meanSpecNeeds_pool[it] = 0
        else
            meanEats_pool[it] = sum(length(poolnet[id].eat) for id in poolnet.spec)/numSpec_pool;     #bit afraid of overflows...
            meanNeeds_pool[it] = sum(length(poolnet[id].need) for id in poolnet.spec)/numSpec_pool;
            meanSpecEats_pool[it] = sum(count(poolnet.hasspec[collect(poolnet[id].eat)]) for id in poolnet.spec)/numSpec_pool;     #bit afraid of overflows...
            meanSpecNeeds_pool[it] = sum(count(poolnet.hasspec[collect(poolnet[id].need)]) for id in poolnet.spec)/numSpec_pool;
        end
        if createLog     #could be significantly optimized by just saving the changes and reconstructing if necessary
            vertsInColony[:,it] = colnet.hasv;
            maxids[it] = poolnet.idmanager.maxid;
            if it >= 9500 && it % 10 == 0
                eatMatrix = ENIgMaGraphs.convertToEatMatrixNonReduced(colnet)
                inds = sort!(vcat(collect(colnet.basalRes),getConnectedSpec(colnet)))

                eatMatrix = eatMatrix[inds,inds]
                R"""
                    library(MASS)  
                    library(NetIndices)
                    rtl<-TrophInd(t($eatMatrix))
                """
                @rget rtl;
                push!(trophLevels, rtl[21:end,:TL] .- 1)

                #oldColonies[it] = deepcopy(colnet)
                #oldPools[it] = deepcopy(poolnet)
            end
        end

        specRich[it] = getNumSpec(colnet);
        rich[it] = length(colnet);
        pool[it] = length(poolnet);

        nPrimExtSpec[it] = length(primextspec)
        nSecExtSpec[it] = length(secextspec)
        nColonizers[it] = length(colonizers)

        #NOTE: standardize mean strength between 0 and 1
    end #end time steps

    if diverse == 1.5
        #finalize phylogenetic tree
        for survivorId in poolnet.spec
            survivingNode = getnode(phyloTree,"$(survivorId)")
            parent = getnode(phyloTree,survivingNode.data["parentName"])
            createbranch!(phyloTree,parent,survivingNode,t - parent.data["timestamp"])
            survivingNode.data["extinct"] = false
            survivingNode.data["timestamp"] = t
        end
    end

    if createLog
        return ENIgMaSimulationData_v3(
            poolnet,
            colnet,
            phyloTree,
            specRich,
            rich,
            pool,
            mstrength,
            clock,
            vertsInColony,
            maxids,
            globextspec,
            meanEats,
            meanNeeds,
            meanEats_pool,
            meanNeeds_pool,
            meanSpecEats,
            meanSpecNeeds,
            meanSpecEats_pool,
            meanSpecNeeds_pool,
            events,
            trophLevels,
            nPrimExtSpec,
            nSecExtSpec,
            nColonizers
        ), 
        ()#oldColonies = oldColonies,oldPools = oldPools)  # for other data that is to be transfered only temporarily
    else
        return ENIgMaSimulationData_v3(
            poolnet,
            colnet,
            phyloTree,
            specRich,
            rich,
            pool,
            mstrength,
            clock,
            BitVector(),
            Int[],
            Dict{Int,Pair{Int,ENIgMaVert}}(),
            meanEats,
            meanNeeds,
            meanEats_pool,
            meanNeeds_pool,
            meanSpecEats,
            meanSpecNeeds,
            meanSpecEats_pool,
            meanSpecNeeds_pool,
            events,
            [],
            nPrimExtSpec,
            nSecExtSpec,
            nColonizers
        ),
        ()#oldColonies = oldColonies,oldPools = oldPools)
    end
end

#keep old versions of simulation Data in order to be able to load them from files
struct ENIgMaSimulationData_v1 <: ENIgMaSimulationData
    poolnet::ENIgMaGraph
    colnet::ENIgMaGraph
    phyloTree::ManyRootTree
    sprich::Vector{Int}
    rich::Vector{Int}
    pool::Vector{Int}
    mstrength::Vector{Float64}
    clock::Vector{Float64}
    CID::BitArray
    maxids::Vector{Int}
    globextspec::Dict{Int,Pair{Int,ENIgMaVert}}
    freqe::Vector{Float64}
    freqn::Vector{Float64}
    freqe_pool::Vector{Float64}
    freqn_pool::Vector{Float64}
    events::Vector{AbstractENIgMaEvent}
end

"""
    Data type that stores the generated data of a simulation of the ENIgMaModel.

# Extended help
    #Fields
    -'poolnet::ENIgMaGraph': The pool network at the end of the simulation.
    -'colnet::ENIgMaGraph': The colony network at the end of the simulation.
    -'phyloTree::ManyRootTree' The phylogenetic tree generated if diversification is on.
    -'specRich::Vector{Int}': The colony's species richness in each itteration.
    -'rich::Vector{Int}': The number of vertices in the colony for each itteration.
    -'pool::Vector{Int}': The number of vertices in the pool for each itteration.
    -'mstrength::Vector{Float64}': The (somewhat normalized) mean strength for each itteration.
    -'clock::Vector{Float64}': The time since the start of the simulation for each itteration.
    -'vertsInColony::BitArray': Stores which vertices are in the colony for each itteration.
        vertsInColony[id,it] is true if vertex with id 'id' is in colony at itteration 'it'.
    -'maxids::Vector{Int}': Saves the maximal vertex id for each itteration.
    -'globextspec::Dict{Int,Pair{Int,ENIgMaVert}}': Stores all (id => vertex) pairs
        of globally extinct species using the itteration they went extinct as keys.
    -'meanEats::Vector{Float64}': The mean out-degree of the eat-subnetwork of the colony.
    -'meanNeeds::Vector{Float64}': The mean out-degree of the need-subnetwork of the colony.
    -'meanEats_pool::Vector{Float64}': The mean out-degree of the eat-subnetwork of the pool.
    -'meanNeeds_pool::Vector{Float64}': The mean out-degree of the need-subnetwork of the pool.
    -'meanSpecEats::Vector{Float64}': The mean out-degree of the species-species eat-subnetwork of the colony.
    -'meanSpecNeeds::Vector{Float64}': The mean out-degree of the species-species need-subnetwork of the colony. 
    -'meanSpecEats_pool::Vector{Float64}': The mean out-degree of the species-species eat-subnetwork of the pool.
    -'meanSpecNeeds_pool::Vector{Float64}': The mean out-degree of the species-species need-subnetwork of the pool. 
    -'events::Vector{AbstractENIgMaEvent}': Stores the event of each itteration.
"""
struct ENIgMaSimulationData_v2 <: ENIgMaSimulationData
    poolnet::ENIgMaGraph
    colnet::ENIgMaGraph
    phyloTree::ManyRootTree
    specRich::Vector{Int}
    rich::Vector{Int}
    pool::Vector{Int}
    mstrength::Vector{Float64}
    clock::Vector{Float64}
    vertsInColony::BitArray
    maxids::Vector{Int}
    globextspec::Dict{Int,Pair{Int,ENIgMaVert}}
    meanEats::Vector{Float64}
    meanNeeds::Vector{Float64}
    meanEats_pool::Vector{Float64}
    meanNeeds_pool::Vector{Float64}
    meanSpecEats::Vector{Float64}
    meanSpecNeeds::Vector{Float64}
    meanSpecEats_pool::Vector{Float64}
    meanSpecNeeds_pool::Vector{Float64}
    events::Vector{AbstractENIgMaEvent}
end