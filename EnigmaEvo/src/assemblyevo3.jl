function assemblyevo(poolnet::ENIgMaGraph, rates, maxits, cm, cn, ce, cf, diverse, restrict_colonization::Bool, createlog = true)
    #create the colony network, use idmanager of pool network (no copy, a reference to the original)
    colnet::ENIgMaGraph = ENIgMaGraph(poolnet.estsize,poolnet.idmanager);
    #add basal resource as its always there
    basalres = ENIgMaVert();
    #add itself as need in order to  make sure it doesnt go extinct (needs one need as object...easier to program that way)
    addn!(basalres,1);
    addmod!( colnet, 1, basalres);
    
    # set up some monitoring buffers
    sprich = Array{Int64}(undef,maxits);
    rich = Array{Int64}(undef,maxits);
    pool = Array{Int64}(undef,maxits);
    mstrength = Array{Float64}(undef,maxits);
    clock = Array{Float64}(undef,maxits);
    events = Array{AbstractENIgMaEvent}(undef,maxits);
    if createlog
        CID = falses(poolnet.estsize,maxits);
        maxids = zeros(Int,maxits);
        globextspec = Dict{Int,Pair{Int,ENIgMaVert}}()    #Stores globally extinct species together with their ids. The keys are the iterations the species went extinct
    end

    #initialize phylogenetic tree
    phyloTree = ManyRootTree()  #create phylo genetic tree
    if diverse == 1  #in non diverse case phylogeny kinda weird
        ##create root nodes
        #createnodes!(phyloTree,Dict{String,Dict{String,Any}}("$(id)_root" => Dict{String,Any}("timestamp" => 0) for id in poolnet.spec))
        superRoot = createnode!(phyloTree,"superRoot", data=Dict{String,Any}("timestamp"=>0.0,"evolution" => ""))
        #create leaf nodes those will be kept and move along with time, see note in evlution passage
        createnodes!(phyloTree,Dict{String,Dict{String,Any}}("$id" => Dict{String,Any}("parentName"=>"superRoot", "version"=>1, "heritage"=>"") for id in poolnet.spec))
    end

    freqe = Array{Float64}(undef,maxits);
    freqn = Array{Float64}(undef,maxits);
    freqe_pool = Array{Float64}(undef,maxits);
    freqn_pool = Array{Float64}(undef,maxits);

    mutstep = Float64[]#zeros(Float64,maxits);
    evolvedstrength = Array{Float64}(undef,0);

    #evolutiontable = [[0 0 0 1 1 1 2 2 2 3 3 3];[1 2 3 0 2 3 0 1 3 0 1 2]];
    #tallytable = [4.1 4.2 5.1 4.3 4.4 5.2 4.5 4.6 5.3 6.1 6.2 6.3];
    interactionTable = Dict{Int, String}(0=>"i",1=>"e",2=>"n",3=>"m")

    # If diversification is turned off, rates.rext -> 0 do that by hand
    if diverse == 0
        rates = (rc = rates.rc,rsecext = rates.rsecext, rprimext = rates.rprimext, reo = rates.reo,revo = rates0.revo,rext = 0.);
    end

    colonizers = Int[]; #later used to store ids of potential colonizers
    secextspec = Int[]; #later used to store ids of spec that could go extinct secondarily
    extobj = Int[];     #later used to store ids of objects that could go extinct

    t=0;
    it = 0;
    while it < maxits

        minstrength = ce*Float64(length(poolnet) - 1);    #assuming ce > cf, eat everything but yourself
        maxstrength = cm*Float64(length(poolnet) - numspec(poolnet) - 1) + cn*Float64(S) - ce*1 - cf*0;  #assuming cm > cn

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
        levo = max(numspec(colnet) - 2,0); #means evolution can only occur if community has more than 2 species

        #COUNT POTENTIAL global extinction events (size of pool)
        lext =  numspec(poolnet) - numspec(colnet) - (colnet.hasv[2] ? 0 : 1) #species 2 can't go exinct

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
            spmutid = rand(collect(setdiff(colnet.spec,2)));    #might be optimized, by finding smarter way, to draw from Set (a#same below for intmutid)
            #OR Mutate realized interactions
            intmutid = rand(collect(setdiff(keys(colnet.vert),[spmutid,2])));    #could possibly be optimized by drawing in loop until valid value is drawn

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
                if intmutid == 1
                    new_int = rand(setdiff(spec_ints,old_int)); #cant make primal resource (sun)
                else
                    new_int = rand(setdiff(object_ints,old_int));
                end
            end

            newId = mutate!(poolnet, colnet, spmutid, intmutid, change_in_int, old_int, new_int, diverse, ce, cn, cm, cf);
            


            if createlog
                if poolnet.estsize > size(CID)[1];      #has the network grown bigger than our buffer?
                    CID = vcat(CID,falses(poolnet.estsize - size(CID)[1],maxits));  #then extend the buffer
                end
            end

            events[it] = MutationEvent(spmutid,intmutid,newId,InteractionType(old_int),InteractionType(new_int), change_in_int);

            #update phylogenic Tree (approach a bit unintuitive.
            # I keep the leaf nodes unattached, just move them forward and let them save their current version
            # and where they are currently attached without actually attaching them)
            if diverse == 1
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
            globextid = rand(collect(setdiff(poolnet.spec,union(colnet.spec,2))));    #could possibly be optimized by finding better way to sample from set
            
            if colnet.hasv[globextid]
                delv!(colnet,globextid)
            end

            if createlog
                globextspec[it] = Pair(globextid,poolnet[globextid]);
            end

            delv!(poolnet,globextid);

            events[it] = GlobalExtinctionEvent(globextid) 

            #update Phylogenetic tree
            if diverse == 1
                nowExtNode = getnode(phyloTree,"$(globextid)")
                parent = getnode(phyloTree,nowExtNode.data["parentName"])
                createbranch!(phyloTree,parent,nowExtNode,t - parent.data["timestamp"])
                nowExtNode.data["extinct"] = true
                nowExtNode.data["timestamp"] = t
            end
        end
        
        #SAVE RESULTS IN BUFFERS
        freqe[it] = (sum([length(colnet[id].eat) for id in colnet.spec]) - length(colnet[1].feed))/max(length(colnet)-1,1);     #bit afraid of overflows...
        freqn[it] = sum([length(setdiff(colnet[id].need,1)) for id in colnet.spec])/max(length(colnet)-1,1);
        freqe_pool[it] = (sum([length(poolnet[id].eat) for id in poolnet.spec]) - length(poolnet[1].feed))/max(length(poolnet)-1,1);     #bit afraid of overflows...
        freqn_pool[it] = sum([length(setdiff(poolnet[id].need,1)) for id in poolnet.spec])/max(length(poolnet)-1,1);

        if createlog     #could be significantly optimized by just saving the changes and reconstructing if necessary
            CID[:,it] = colnet.hasv;
            maxids[it] = poolnet.idmanager.maxid;
        end

        sprich[it] = numspec(colnet);
        rich[it] = length(colnet);
        pool[it] = length(poolnet);
        #NOTE: standardize mean strength between 0 and 1
        mstrength[it] = (mean( v.strength for (id,v) in poolnet if poolnet.hasspec[id] ) - minstrength)/(maxstrength - minstrength); #make strengths positive with a minimum of 1
    end #end time steps

    if diverse == 1
        #finalize phylogenetic tree
        for survivorId in poolnet.spec
            survivingNode = getnode(phyloTree,"$(survivorId)")
            parent = getnode(phyloTree,survivingNode.data["parentName"])
            createbranch!(phyloTree,parent,survivingNode,t - parent.data["timestamp"])
            survivingNode.data["extinct"] = false
            survivingNode.data["timestamp"] = t
        end
    end

    if createlog
        return(
            poolnet,
            colnet,
            phyloTree,
            sprich,
            rich,
            pool,
            mstrength,
            evolvedstrength,
            clock,
            CID,
            maxids,
            globextspec,
            mutstep,
            freqe,
            freqn,
            freqe_pool,
            freqn_pool,
            events
        )
    else
        return(
            poolnet,
            colnet,
            phyloTree,
            sprich,
            rich,
            pool,
            mstrength,
            evolvedstrength,
            clock,
            BitVector(),
            Int[],
            Dict{Int,Pair{Int,ENIgMaVert}}(),
            mutstep,
            freqe,
            freqn,
            freqe_pool,
            freqn_pool,
            events
        )
    end
end