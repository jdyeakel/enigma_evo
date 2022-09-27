function assemblyevo(poolnet::ENIgMaGraph, rates0, maxits, cm, cn, ce, cf, diverse)
    # Total size of the species + objects
    N = length(poolnet);
    S = numspec(poolnet);
    #Number of objects
    O = N - S - 1;

    #create the colony network with an estimated max size (sizehint) dependent on wheter or not diversification is activated
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
    CID = falses(poolnet.estsize,maxits);

    freqe = Array{Float64}(undef,maxits);
    freqn = Array{Float64}(undef,maxits);

    mutstep = zeros(Float64,maxits);
    evolvedstrength = Array{Float64}(undef,0);

    evolutiontable = [[0 0 0 1 1 1 2 2 2 3 3 3];[1 2 3 0 2 3 0 1 3 0 1 2]];
    tallytable = [4.1 4.2 5.1 4.3 4.4 5.2 4.5 4.6 5.3 6.1 6.2 6.3];

    # If diversification is turned off, rates.rext -> 0
    if diverse == 1
        rates = (rc = rates0.rc,re = rates0.re,reo = rates0.reo,revo = rates0.revo,rext = rates0.rext);
    else
        rates = (rc = rates0.rc,re = rates0.re,reo = rates0.reo,revo = rates0.revo,rext = 0.);
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
        getpotcolonizers!(poolnet,colnet,colonizers);   #saves pot colonizers ids in colonizers

        #COUNT SPECIES THAT COULD POTENTIALLY GO EXTINCT
        
        #1) SECONDARY EXTINCTIONS By not fulfilling Eat/Need thresholds
        getsecext!(colnet,secextspec,extobj);

        #2) PRIMARY EXTINCTIONS by competition
        primextspec = getprimext(colnet,ce,cn,cm,cf);

        #combine both extinction types
        specxt = union(primextspec,secextspec)
        lspext = length(specxt)

        #COUNT POTENTIAL evolutionary events (size of community)
        levo = maximum([numspec(colnet) - 2,0]); #means evolution can only occur if community has more than >2 species

        #COUNT POTENTIAL global extinction events (size of pool)
        lext = max(numspec(colnet)-1,0); #-1 to account for species 1, which is restricted


        # Calculate the full Rate
        Rate = rates.rc*length(colonizers) + rates.re*lspext + rates.reo*length(extobj) + rates.revo*levo + rates.rext*lext;

        # Calculate event probabilities
        probc = (rates.rc*length(colonizers))/Rate;
        probe = (rates.re*lspext)/Rate;
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
        if re < probc #(lcol/levents)

            #COLONIZATION FUNCTION
            idc = rand(colonizers);
            colonize!(poolnet,colnet,idc)

            #tally
            tally = 0;

        #DRAW EXTINCTION
        elseif re < (probc + probe)

            #SPECIES EXTINCTION FUNCTION
            #select species to go extinct
            sp_bye = rand(specxt);

            #if sp_bye == 341
            #    println("debugging")
            #end

            delv!( colnet, sp_bye )

            #tally
            if in(sp_bye,primextspec) 
                #These are extinctions from competitive exclusion
                tally = 1;
            end
            if in(sp_bye,secextspec)
                #These are secondary extinctions
                tally = 2;
            end

        #DRAW OBJECT EXTINCTION
        elseif re < (probc + probe + probeo) #((lcol + lspext)/levents) && re < ((lcol + lspext + lobext)/levents)

            #OBJECT EXTINCTION FUNCTION
            delv!( colnet, rand(extobj) )

            #tally
            tally = 3;
        
        #DRAW EVOLUTION
        elseif re < (probc + probe + probeo + probevo)

            # EVOLUTION OF SPECIES IN THE COMMUNITY
            # SPECIES MUTATION
            # select FROM community
            spmutid = rand(collect(setdiff(colnet.spec,2)));    #might be optimized, by finding smarter way, to draw from Set (a#same below for intmutid)
            #OR Mutate realized interactions
            intmutid = rand(collect(setdiff(keys(colnet.vert),spmutid)));    #could possibly be optimized by drawing in loop until valid value is drawn

            spints = [0,1,2];
            obints = [0,1,2,3];
            tally = mutate!(poolnet, colnet, spmutid, intmutid, diverse, spints, obints, tallytable, evolutiontable,ce,cn,cm,cf);

        # Global EXTINCTION
        else
            
            #Choose species subject to extinction
            #Do not choose species 2, which is protected
            globextid = sample(collect(setdiff(poolnet.spec,2)));    #could possibly be optimized by finding better way to sample from set
            
            if colnet.hasv[globextid]
                delv!(colnet,globextid)
            end
            
            #update counters
            if poolnet.hasspec[globextid]   #is species?
                S -= 1;
            else
                O -= 1;
            end

            delv!(poolnet,globextid);
            N -= 1;
        end
        
        #SAVE RESULTS IN BUFFERS
        freqe[it] = sum([length(vert.eat) for (_,vert) in colnet])/length(colnet.vert);     #bit afraid of overflows...
        freqn[it] = sum([length(vert.need) for (_,vert) in colnet])/length(colnet.vert);

        if diverse == 0     #could be optimized (improved, not performance wise) by including the diversification case... Would just have to add CID to a list of objects that might eventually need to be extended
            CID[:,it] = colnet.hasspec;
        end

        sprich[it] = numspec(colnet);
        rich[it] = length(colnet);
        pool[it] = N;
        #NOTE: standardize mean strength between 0 and 1
        mstrength[it] = (mean( v.strength for (id,v) in poolnet if poolnet.hasspec[id] ) - minstrength)/(maxstrength - minstrength); #make strengths positive with a minimum of 1
        events[it] = tally;
    end #end time steps


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
        #intm_evo,
        mutstep,
        freqe,
        freqn,
        events
    )

end
