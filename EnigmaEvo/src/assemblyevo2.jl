function assemblyevo(S,intm,eb,nb,nb0,mb,e_t,n_t,maxits,probmut,cm,cn,ce,cpred)

    # S = length(spv) + 1;
    # Total size of the species + objects
    N = size(intm)[1];
    # A list of the species ids
    spv = collect(Int64,2:1:S);
    spobv = collect(Int64,1:1:N);

    # Make a copy of the original interaction matrix
    intm_orig = copy(intm);

    # MaxN = convert(Int64,floor(S + S*lambda));
    cid = Array{Int64}(undef,0);

    sprich = Array{Int64}(undef,maxits);
    rich = Array{Int64}(undef,maxits);
    mstrength = Array{Float64}(undef,maxits);
    clock = Array{Float64}(undef,maxits);
    events = Array{Float64}(undef,maxits);
    CID = falses(N,maxits);

    freqe = Array{Float64}(undef,maxits);
    freqn = Array{Float64}(undef,maxits);

    mutstep = zeros(Float64,maxits);
    evolvedstrength = Array{Float64}(undef,0);
    #NOTE strength matrix can't be built a-priori!
    # #Build the strength matrix apriori
    # strength = vec(pi*sum(nb0,dims=2)) .- vec(sqrt(2)*sum(eb,dims=2)) .- vec(sum(eb,dims=1));
    # smatrix = Array{Float64}(copy(eb));
    # for i=1:size(eb)[1]
    #     smatrix[i,:] .= eb[i,:] * strength[i];
    # end
    # # smatrix[findall(iszero,eb)] = NaN;
    #
    minstrength = 0 -ce*Float64(N) - cpred*Float64(S);
    maxstrength = cn*Float64(N) - ce*1 -cpred*0;

    evolutiontable = [[0 0 0 1 1 1 2 2 2 3 3 3];[1 2 3 0 2 3 0 1 3 0 1 2]];
    tallytable = [4.1 4.2 5.1 4.3 4.4 5.2 4.5 4.6 5.3 6.1 6.2 6.3];


    t=0;
    it = 0;
    while it < maxits

        #does sorting CID make a difference?
        sort!(cid);

        cid_old = copy(cid);
        #Which are species?
        spcid = intersect(spv,cid);
        #Which are objects?
        ocid = setdiff(cid,spcid);


        #COUNT POTENTIAL COLONIZERS
        trophiclinked = setdiff(spobv[vec(sum(eb[:,[1;cid]],dims=2) .> 0)],cid);
        #For each trophiclinked, count number of assimilate and need interactions in system
        #Determine in the proportion that already exists is >= the threshold
        e_fill = ((sum(eb[trophiclinked,[1;cid]],dims=2)./sum(eb[trophiclinked,:],dims=2)) .>= e_t)[:,1];
        prop_n = sum(nb0[trophiclinked,[1;cid]],dims=2)./sum(nb0[trophiclinked,:],dims=2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] .= 1;
        n_fill = (prop_n .>= n_t)[:,1];
        #Build a list of those that pass each individual test
        e_pass = trophiclinked[e_fill];
        n_pass = trophiclinked[n_fill];
        #Build a list of those that pass both tests (potential colonizers)
        col = intersect(e_pass,n_pass);
        #NOTE: alternative: Free colonization with one eat.
        # col = copy(e_pass);
        #Count the number that pass
        lcol = length(col);


        #COUNT POTENTIAL EXTINCT SPECIES
        #1) SECONDARY EXTINCTIONS By not fulfilling Eat/Need thresholds
        #Re-calculate e_t and n_t (> e_t; >= n_t)
        secext = secexteval(spcid,cid,eb,nb0,e_t,n_t);

        if length(spcid) > 0

            # Build the strength matrix at each community state
            # Strength values change over time so need to ve updated
            # Only record strength values of species (hence the [1:length(spcid)])
            #NOTE: Needs won't change; Eats is based on POTENTIAL niche; Vuln changes per timestep
            
            #NOTE: trying a new strength function here
            # strength = vec(cn*sum(nb0[spcid,cid],dims=2)) .- vec(ce*sum(eb[spcid,:],dims=2)) .- (vec(cpred*sum(eb[spcid,cid],dims=1))[1:length(spcid)]);
            
            strength = Array{Float64}(undef,length(spcid));
            cmatrix = Array{Float64}(undef,length(spcid),length(cid));
            for i=1:length(spcid)
                strength[i] = strengthcalc(nb0,eb,cid,spcid[i],cm,cn,ce,cpred);
                cmatrix[i,:] .= eb[spcid[i],cid] * strength[i];
            end

            # cmatrix = eb[spcid,cid] .* reshape(repeat(strength,outer=length(cid)),length(spcid),length(cid));

            # cmatrix = (eb[spcid,cid]' * (Matrix{Float64}(I,length(strength),length(strength)) .* strength))'

            #'zero' entrees need to be lower than any possible strength
            #So they will be effectively ignored
            #-sqrt(2)*S-S is the theoretical min strength
            cmatrix[cmatrix.==0] .= minstrength;

            #Maximum strength values for each resource utilized
            #Competitor has to match this value for at least one of its foods
            cmax = findmax(cmatrix,dims=1)[1];

            prext_comp = trues(length(spcid));
            for i=1:length(spcid)

                #Don't count sun
                #catalogue prey for all species/objects in the system (excluding sun)
                ieats = Array{Bool}(eb[spcid[i],cid]);
                #If you have >= the max strength for any of those prey, you stay
                #This means that a pure primary producer is not evaluated

                # any(ieats) is false if you are a pure primary producer, meaning you can't go extinct from competition in this model
                # (any(strength[i] .>= cmax[ieats])==false) is false if you are the strongest (or equally strongest) competitor among consumers that share your food
                prext_comp[i] = any(ieats)*(any(strength[i] .>= cmax[ieats])==false);

            end
            primext = spcid[prext_comp];
        else
            strength = minstrength;
            primext = Array{Int64}(undef,0);
        end

        # primext = copy(spext2);
        # secext = copy(spext1);
        #
        # spext2 = Array{Int64}(0);
        spext = unique([secext;primext]);
        lspext = length(spext);



        #COUNT POTENTIAL EXTINCT OBJECTS
        obext = ocid[findall(iszero,vec(sum(mb[spcid,ocid],dims=1)))];
        lobext = length(obext);


        # lmutations = length(cid);

        lecoevents = sum([lcol;lspext;lobext]);
        if length(spcid) < 2
            probevo = 0.;
        else
            probevo = copy(probmut);
        end
        #Total number of events include ECO events + EVO event
        levents = Int64(floor((lecoevents)/(1-probevo))); 
       
        # #mutation probability is set (and not dependent on state)
        # lmutations = Int64(floor(probmut*levents));
        # #No mutations if there are no interactions in the system
        # if length(spcid) < 2
        #     lmutations = 0;
        # end
        # #Redefine levents to include lmutations
        # levents = sum([lcol;lspext;lobext;lmutations]);

        dt = 1/levents;

        t += dt;
        it += 1;
        # push!(clock,t);
        clock[it] = t;

        #Choose a random event
        re = rand();
        tally = -10;

        #DRAW COLONIZATION
        if re < (lcol/levents)

            #COLONIZATION FUNCTION
            c_sp = rand(col);
            #select made objects as well
            c_ob_full = findall(isodd,mb[c_sp,:]);
            #Only bring in objects that aren't already in the community
            c_ob = setdiff(c_ob_full,cid_old);
            colonize = [c_sp;c_ob];
            #Update community
            cid = [cid_old;colonize];
            #Update species
            spcid = [spcid;c_sp];
            #tally
            tally = 0;
        end

        #DRAW EXTINCTION
        if re > (lcol/levents) && re < ((lcol + lspext)/levents)

            #SPECIES EXTINCTION FUNCTION
            #select species to go extinct
            sp_bye = rand(spext,1);
            cid = setdiff(cid_old,sp_bye);
            spcid = setdiff(spcid,sp_bye);
            #tally
            if in(sp_bye[1],primext)
                #These are extinctions from competitive exclusion
                tally = 1;
            end
            if in(sp_bye[1],secext)
                #These are secondary extinctions
                tally = 2;
            end
        end

        #DRAW OBJECT EXTINCTION
        if re > ((lcol + lspext)/levents) && re < ((lcol + lspext + lobext)/levents)

            #OBJECT EXTINCTION FUNCTION
            ob_bye = rand(obext,1);
            cid = setdiff(cid_old,ob_bye);
            #tally
            tally = 3;
        end
        
        #DRAW EVOLUTION
        if re > ((lcol + lspext + lobext)/levents)

            re2 = 0.; #rand();
            lcom = length(cid);
            lpool = N - lcom;
            levents2 = lcom + lpool;

            # EVOLUTION OF SPECIES IN THE COMMUNITY
            if re2 < lcom/levents2

                #SPECIES MUTATION
                #select FROM COMMUNITY
                spmut = rand(spcid);
                #OR Mutate realized interactions
                intmut = rand(setdiff(cid,spmut));

                spints = [0,1,2];
                obints = [0,1,2,3];
                intm_mut, ebmut, nbmut, nb0mut, mbmut, tally = mutation(spmut,intmut,spints,obints,intm,spv,tallytable,evolutiontable);

                #Does the mutant outcompete the parent?
                strengthmut = strengthcalc(nb0mut,ebmut,cid,spmut,cm,cn,ce,cpred);
                #Evaluate if spmut will go secondarily extinct (through disconnection)
                # secextmut = secexteval(zeros(Int64,1).+spmut,cid,ebmut,nb0mut,e_t,n_t);
                #If it does probability is 1
                # probsecextmut = (length(secextmut) > 0)*1.;
                # probsecextmut = 0.
                
                spmutpos = findall(x->x==spmut,spcid)[1];
                muteats = ebmut[spmut,cid] .== 1;
                spmutstrength = muteats * strengthmut;
                spmutstrength[spmutstrength.==0] .= minstrength;
                # Original species competition strength
                # spstrength = cmatrix[spmutpos,:];
                
                if any(spmutstrength[muteats] .> cmax[muteats]) #&& probsecextmut != 1.
                    #accept mutation
                    intm = copy(intm_mut);
                    eb = copy(ebmut);
                    nb = copy(nbmut);
                    nb0 = copy(nb0mut);
                    mb = copy(mbmut)
                    # i_b = copy(i_bmut);

                    mutstep[it] = 1.0
                    push!(evolvedstrength, strengthmut - strength[spmutpos]);

                end
            else
                
                # RANDOM EVOLUTION OF SPECIES IN THE POOL
                #species pool
                sppool = setdiff(spv,cid);
                #Object pool
                obpool = setdiff(spobv,cid);
                spmut = rand(sppool);
                intmut = rand(setdiff([sppool;obpool],spmut));

                # Random mutation
                spints = [0,1,2];
                obints = [0,1,2,3];
                # spints = vec(((eb .* 1) .+ (nb0 .* 2))[2:S,1:S]);
                # obints = vec(((eb .* 1) .+ (nb0 .* 2))[2:S,(S+1):N]);
                intm_mut, ebmut, nbmut, nb0mut, mbmut, tally = mutation(spmut,intmut,spints,obints,intm,spv,tallytable,evolutiontable);

                # NOTE: this section doesn't work yet (maybe not necessary)
                # Mutate back to original
                # diffsp = findall(isodd,vec(sum(intm[[sppool;obpool],:] .!= intm_orig[[sppool;obpool],:],dims=2)));
                # evolvedinteractions = findall(isodd,);
                # rand(collect(1:length(evolvedinteractions)));
                
                #accept mutation if spmut still has something it consumes (don't allow it to de-evolve interactions completely)
                if sum(eb_mut[spmut,:]) >= 1
                    intm = copy(intm_mut);
                    eb = copy(ebmut);
                    nb = copy(nbmut);
                    nb0 = copy(nb0mut);
                    mb = copy(mbmut)
                    # i_b = copy(i_bmut);
                    mutstep[it] = 2.0;
                end

            end

        end

        freqe[it] = sum(eb[spcid,cid])/(length(cid));
        freqn[it] = sum(nb0[spcid,cid])/(length(cid));

        #NOTE - updating CID....
        CID[cid,it] .= true;
        sprich[it] = length(spcid);
        rich[it] = length(cid);
        #NOTE: standardize mean strength between 0 and 1
        mstrength[it] = (mean(strength) - minstrength)/(maxstrength - minstrength); #make strengths positive with a minimum of 1
        events[it] = tally;
        # push!(sprich,length(spcid));
        # push!(rich,length(cid));
        # push!(events,tally);
    end #end time steps

    intm_evo = copy(intm);



    return(
    sprich,
    rich,
    mstrength,
    evolvedstrength,
    clock,
    CID,
    intm_evo,
    mutstep,
    freqe,
    freqn,
    events
    )

end
