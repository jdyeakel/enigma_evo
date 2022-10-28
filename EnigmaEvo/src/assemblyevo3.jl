function assemblyevo(poolnet::ENIgMaGraph, rates, maxits, cm, cn, ce, cf, diverse, restrict_colonization::Bool, createlog = false)
    # Total size of the species + objects
    N = length(poolnet);
    S = numspec(poolnet);
    #Number of objects
    O = N - S - 1;

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
    events = Array{Float64}(undef,maxits);
    if createlog
        CID = falses(poolnet.estsize,maxits);
        maxids = zeros(Int,maxits);
        globextspec = Dict{Int,Pair{Int,ENIgMaVert}}()    #Stores globally extinct species together with their ids. The keys are the iterations the species went extinct
    end

    freqe = Array{Float64}(undef,maxits);
    freqn = Array{Float64}(undef,maxits);

    mutstep = zeros(Float64,maxits);
    evolvedstrength = Array{Float64}(undef,0);

    evolutiontable = [[0 0 0 1 1 1 2 2 2 3 3 3];[1 2 3 0 2 3 0 1 3 0 1 2]];
    tallytable = [4.1 4.2 5.1 4.3 4.4 5.2 4.5 4.6 5.3 6.1 6.2 6.3];

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

        minstrength = 0 + 0 - ce*Float64(N - 1) - cf*Float64(O);    #assuming ce > cf
        maxstrength = cm*Float64(O) + cn*Float64(S) - ce*1 - cf*0;  #assuming cm > cn

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
        tally = -10;

        #DRAW COLONIZATION
        if re <= probc #(lcol/levents)

            #COLONIZATION FUNCTION
            idc = rand(colonizers);
            colonize!(poolnet,colnet,idc)

            #tally
            tally = 0;

        #DRAW primary EXTINCTION
        elseif re <= (probc + probprimext)
            #select species to go extinct
            sp_bye = rand(primextspec);
            delv!( colnet, sp_bye )

            #These are extinctions from competitive exclusion
            tally = 1;

        #DRAW SECONDARY EXTINCTION
        elseif re <= (probc + probprimext + probsecext)

            #select species to go extinct
            sp_bye = rand(secextspec);
            delv!( colnet, sp_bye )

            #These are extinctions from competitive exclusion
            tally = 2;

        #DRAW OBJECT EXTINCTION
        elseif re <= (probc + probprimext + probsecext + probeo) #((lcol + lspext)/levents) && re < ((lcol + lspext + lobext)/levents)

            #OBJECT EXTINCTION FUNCTION
            delv!( colnet, rand(extobj) )

            #tally
            tally = 3;
        
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

            tally = mutate!(poolnet, colnet, spmutid, intmutid, change_in_int, old_int, new_int, diverse, tallytable, evolutiontable, ce, cn, cm, cf);
            
            if createlog
                if poolnet.estsize > size(CID)[1];      #has the network grown bigger than our buffer?
                    CID = vcat(CID,falses(poolnet.estsize - size(CID)[1],maxits));  #then extend the buffer
                end
            end

        # Global EXTINCTION
        elseif re <= 1
            #Choose species subject to extinction
            #Do not choose species 2, which is protected
            globextid = rand(collect(setdiff(poolnet.spec,union(colnet.spec,2))));    #could possibly be optimized by finding better way to sample from set
            
            if colnet.hasv[globextid]
                delv!(colnet,globextid)
            end
            
            #update counters
            if poolnet.hasspec[globextid]   #is species?
                S -= 1;
            else
                O -= 1;
            end

            if createlog
                globextspec[it] = Pair(globextid,poolnet[globextid]);
            end

            delv!(poolnet,globextid);
            N -= 1;
        end
        
        #SAVE RESULTS IN BUFFERS
        freqe[it] = (sum([length(colnet[id].eat) for id in colnet.spec]) - length(colnet[1].feed))/max(length(colnet)-1,1);     #bit afraid of overflows...
        freqn[it] = sum([length(setdiff(colnet[id].need,1)) for id in colnet.spec])/max(length(colnet)-1,1);

        if createlog     #could be significantly optimized by just saving the changes and reconstructing if necessary
            CID[:,it] = colnet.hasv;
            maxids[it] = poolnet.idmanager.maxid;
        end

        sprich[it] = numspec(colnet);
        rich[it] = length(colnet);
        pool[it] = N;
        #NOTE: standardize mean strength between 0 and 1
        mstrength[it] = (mean( v.strength for (id,v) in poolnet if poolnet.hasspec[id] ) - minstrength)/(maxstrength - minstrength); #make strengths positive with a minimum of 1
        events[it] = tally;
    end #end time steps

    if createlog
        return(
            poolnet,
            colnet,
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
            events
        )
    else
        return(
            poolnet,
            colnet,
            sprich,
            rich,
            pool,
            mstrength,
            evolvedstrength,
            BitVector(),
            Int[],
            Dict{Int,Pair{Int,ENIgMaVert}}(),
            clock,
            mutstep,
            freqe,
            freqn,
            events
        )
    end
end