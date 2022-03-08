function assemblyevo(rates0,S,intm,eb,nb,nb0,mb,e_t,n_t,maxits,cm,cn,ce,cpred,diverse)

    # S = length(spv) + 1;
   
    # A list of the species ids

    #These will be updated if diverse = 1
    # Total size of the species + objects
    N = size(intm)[1];
    #Number of objects
    O = N-S;
    spv = collect(Int64,2:1:S);
    spobv = collect(Int64,2:1:N); #NOTE: fixed this to start at 2 on 3/8/2022
    espv = Array{Int64}(undef,0);

    # Make a copy of the original interaction matrix
    # intm_orig = copy(intm);

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
    # rates = (rc = 1., re = 1., reo = 1., revo = 1., rext = 1.);

   

    evolutiontable = [[0 0 0 1 1 1 2 2 2 3 3 3];[1 2 3 0 2 3 0 1 3 0 1 2]];
    tallytable = [4.1 4.2 5.1 4.3 4.4 5.2 4.5 4.6 5.3 6.1 6.2 6.3];


    # If diversification is turned off, rates.rext -> 0
    if diverse == 1
        rates = (rc = rates0.rc,re = rates0.re,reo = rates0.reo,revo = rates0.revo,rext = rates0.rext);
    else
        rates = (rc = rates0.rc,re = rates0.re,reo = rates0.reo,revo = rates0.revo,rext = 0.);
    end

    t=0;
    it = 0;
    while it < maxits

        minstrength = 0 + 0 -ce*Float64(N) - cpred*Float64(S);
        maxstrength = cm*Float64(S) + cn*Float64(N) - ce*1 -cpred*0;

        #does sorting CID make a difference?
        sort!(cid);

        cid_old = copy(cid);
        #Which are species?
        spcid = intersect([spv;espv],cid);
        #Which are objects?
        ocid = setdiff(cid,spcid);


        #COUNT POTENTIAL COLONIZERS
        trophiclinked = setdiff(spobv[vec(sum(eb[spobv,[1;cid]],dims=2) .> 0)],cid);
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
                strength[i] = strengthcalc(nb0,eb,mb,cid,spcid[i],cm,cn,ce,cpred);
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

        # lecoevents = sum([lcol;lspext;lobext]);
        # if length(spcid) < 2
        #     probevo = 0.;
        # else
        #     probevo = copy(probmut);
        # end
        # #Total number of events include ECO events + EVO event
        # levents = Int64(floor((lecoevents)/(1-probevo))); 

        #COUNT POTENTIAL evolutionary events (size of community)
        levo = maximum([length(spcid) - 2,0]); #means evolution can only occur if community has more than >2 species


        #COUNT POTENTIAL global extinction events (size of pool)
        lext = (N - O) - 2; #-2 to account for basal resource + species 2, which is restricted


        # Calculate the full Rate
        Rate = rates.rc*lcol + rates.re*lspext + rates.reo*lobext + rates.revo*levo + rates.rext*lext;

        # Calculate event probabilities
        probc = (rates.rc*lcol)/Rate;
        probe = (rates.re*lspext)/Rate;
        probeo = (rates.reo*lobext)/Rate;
        probevo = (rates.revo*levo)/Rate;
        probext = (rates.rext*lext)/Rate;
        # sum([probc,probe,probeo,probevo,probext])
       
        # #mutation probability is set (and not dependent on state)
        # lmutations = Int64(floor(probmut*levents));
        # #No mutations if there are no interactions in the system
        # if length(spcid) < 2
        #     lmutations = 0;
        # end
        # #Redefine levents to include lmutations
        # levents = sum([lcol;lspext;lobext;lmutations]);

        # dt = 1/levents;
        dt = 1/Rate;
        t += dt;
        it += 1;
        # push!(clock,t);
        clock[it] = t;

        #Choose a random event
        re = rand();
        tally = -10;

        #DRAW COLONIZATION
        if re < probc #(lcol/levents)

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
        if re > probc && re < (probc + probe)#(lcol/levents) && re < ((lcol + lspext)/levents)

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
        if re > (probc + probe) && re < (probc + probe + probeo) #((lcol + lspext)/levents) && re < ((lcol + lspext + lobext)/levents)

            #OBJECT EXTINCTION FUNCTION
            ob_bye = rand(obext,1);
            cid = setdiff(cid_old,ob_bye);
            #tally
            tally = 3;
        end
        
        #DRAW EVOLUTION
        if re > (probc + probe + probeo) && re < (probc + probe + probeo + probevo)#((lcol + lspext + lobext)/levents)

            # re2 = 0.; #rand(); #NOTE: turned off if re2 = 0.
            # lcom = length(cid);
            # lpool = N - lcom;
            # levents2 = lcom + lpool;

            # EVOLUTION OF SPECIES IN THE COMMUNITY
            # if re2 < lcom/levents2

            #SPECIES MUTATION
            #select FROM COMMUNITY
            spmut = rand(spcid);
            #OR Mutate realized interactions
            intmut = rand(setdiff(cid,spmut));

            spints = [0,1,2];
            obints = [0,1,2,3];
            intm_mut, ebmut, nbmut, nb0mut, mbmut, tally = mutation(spmut,intmut,spints,obints,intm,spv,tallytable,evolutiontable);

            #Does the mutant outcompete the parent?
            strengthmut = strengthcalc(nb0mut,ebmut,mbmut,cid,spmut,cm,cn,ce,cpred);
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
            
            if diverse == 0 && any(spmutstrength[muteats] .> cmax[muteats]) #&& probsecextmut != 1.
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
            if diverse == 1
                #Diversify the interaction matrix
                intm_div = Array{Int64}(undef,N+1,N+1);
                intm_div[1:N,1:N] = intm;
                intm_div[N+1,1:N] = intm_mut[spmut,1:N];
                intm_div[1:N,N+1] = intm_mut[1:N,spmut];
                intm_div[N+1,N+1] = 2;

                #How does new species interact with its parent taxon?
                intm_div[N+1,spmut] = 0;
                intm_div[spmut,N+1] = 0;
                

                intm = copy(intm_div);
                eb,nb,nb0,mb = intbool(intm_div);


                push!(espv,N+1);
                push!(spobv,N+1);
                
                S += 1; #species + evolved species
                N += 1; #species + objects + evolved species
                

                #Add that species to the community
                #Note: does not seem to be needed to enable growth
            end
            # else
                
            #     # RANDOM EVOLUTION OF SPECIES IN THE POOL
            #     #species pool
            #     sppool = setdiff(spv,cid);
            #     #Object pool
            #     obpool = setdiff(spobv,cid);
            #     spmut = rand(sppool);
            #     intmut = rand(setdiff([sppool;obpool],spmut));

            #     # Random mutation
            #     spints = [0,1,2];
            #     obints = [0,1,2,3];
            #     # spints = vec(((eb .* 1) .+ (nb0 .* 2))[2:S,1:S]);
            #     # obints = vec(((eb .* 1) .+ (nb0 .* 2))[2:S,(S+1):N]);
            #     intm_mut, ebmut, nbmut, nb0mut, mbmut, tally = mutation(spmut,intmut,spints,obints,intm,spv,tallytable,evolutiontable);

            #     # NOTE: this section doesn't work yet (maybe not necessary)
            #     # Mutate back to original
            #     # diffsp = findall(isodd,vec(sum(intm[[sppool;obpool],:] .!= intm_orig[[sppool;obpool],:],dims=2)));
            #     # evolvedinteractions = findall(isodd,);
            #     # rand(collect(1:length(evolvedinteractions)));
                
            #     #accept mutation if spmut still has something it consumes (don't allow it to de-evolve interactions completely)
            #     if sum(eb_mut[spmut,:]) >= 1
            #         intm = copy(intm_mut);
            #         eb = copy(ebmut);
            #         nb = copy(nbmut);
            #         nb0 = copy(nb0mut);
            #         mb = copy(mbmut)
            #         # i_b = copy(i_bmut);
            #         mutstep[it] = 2.0;
            #     end

            # end
        end
        # Global EXTINCTION
        if re > (probc + probe + probeo + probevo)
            
            #Choose species subject to extinction
            #Do not choose species 2, which is protected
            sp_globalext = rand(setdiff([spv;espv],2),1)[1];
            
            #1) First remove from the local community
            #If its not in the local community, this won't do anything
            #setdiff! appears to be faster than renaming
            setdiff!(cid,sp_globalext);
            setdiff!(spcid,sp_globalext);
            
            #2) update intm
            intm_ext = Array{Int64}(undef,N-1,N-1);
            intm_ext[1:(sp_globalext-1),1:(sp_globalext-1)] = intm[1:(sp_globalext-1),1:(sp_globalext-1)];
            intm_ext[1:(sp_globalext-1),sp_globalext:(N-1)] = intm[1:(sp_globalext-1),(sp_globalext+1):N];
            intm_ext[sp_globalext:(N-1),1:(sp_globalext-1)] = intm[(sp_globalext+1):N,1:(sp_globalext-1)];
            intm_ext[sp_globalext:(N-1),sp_globalext:(N-1)] = intm[(sp_globalext+1):N,(sp_globalext+1):N];

            intm = copy(intm_ext);
            eb,nb,nb0,mb = intbool(intm_ext);

            #3) update labels
            #If an original species have to change S and S-dependencies
            if sp_globalext <= S
                S -= 1;
                #reset original species id vector
                spv = collect(Int64,2:1:S);
            end
            
            #If original or evolved, have to change N and N-dependencies
            N -= 1;
            #NOTE: number of objects O always stays the same (for now)

            #this vector is the same size, but ids reduced by 1
            espv = collect(Int64,(S+O+1):1:N);
            spobv = collect(Int64,2:1:N);

            #4) Rename species in the community according to new positions in intm matrix
            cid[findall(x->x>sp_globalext,cid)] .-= 1;
            spcid[findall(x->x>sp_globalext,spcid)] .-= 1;

            # #Is it an original species?
            # if sp_globalext <= S
                
            #     #1) update intm
            #     intm_ext = Array{Int64}(undef,N-1,N-1);
            #     intm_ext[1:(sp_globalext-1),1:(sp_globalext-1)] = intm[1:(sp_globalext-1),1:(sp_globalext-1)];
            #     if sp_globalext < N
            #         intm_ext[(sp_globalext+1):N,(sp_globalext+1):N] = intm_ext[(sp_globalext+1):N,(sp_globalext+1):N];
            #     end
            #     intm = copy(intm_ext);
            #     eb,nb,nb0,mb = intbool(intm_ext);

            #     #2) update labels
            #     S -= 1;
            #     N -= 1;
            #     #reset original species id vector
            #     spv = collect(Int64,2:1:S);
            #     #this vector is the same size, but ids reduced by 1
            #     espv .-= 1; 
            #     spobv = collect(Int64,2:1:N);

            #     #Rename species in the community according to new positions in intm matrix
            #     cid[findall(x->x>sp_globalext,cid)] .-= 1;
            #     spcid[findall(x->x>sp_globalext,cid)] .-= 1;

            # #Or is it an evolved species?
            # else
            #     #1) update intm
            #     intm_ext = Array{Int64}(undef,N-1,N-1);
            #     intm_ext[1:(sp_globalext-1),1:(sp_globalext-1)] = intm[1:(sp_globalext-1),1:(sp_globalext-1)];
            #     if sp_globalext < N
            #         intm_ext[(sp_globalext+1):N,(sp_globalext+1):N] = intm_ext[(sp_globalext+1):N,(sp_globalext+1):N];
            #     end
            #     intm = copy(intm_ext);
            #     eb,nb,nb0,mb = intbool(intm_ext);

            #     #2) update labels
            #     N -= 1;
            #     espv = collect(Int64,(S+O+1):1:N);
            #     spobv = collect(Int64,2:1:N);

            #     #Rename species in the community according to new positions in intm matrix
            #     cid[findall(x->x>sp_globalext,cid)] .-= 1;
            #     spcid[findall(x->x>sp_globalext,cid)] .-= 1;

            # end
            
        end
        

        freqe[it] = sum(eb[spcid,cid])/(length(cid));
        freqn[it] = sum(nb0[spcid,cid])/(length(cid));

        #NOTE: - updating CID only without diversification dynamics
        #Because the size of CID will change dynamically w/ diversification/global extinction
        if diverse == 0
            CID[cid,it] .= true;
        end
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
