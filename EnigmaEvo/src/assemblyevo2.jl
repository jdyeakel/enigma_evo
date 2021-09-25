function assemblyevo(intm,eb,nb,i_b,m_b,nb0,sp_v,int_id,lambda,
    e_t,n_t,maxits,probmut,cn,ce,cp)

    S = length(sp_v) + 1;
    N = size(intm)[1];
    MaxN = convert(Int64,floor(S + S*lambda));
    cid = Array{Int64}(undef,0);
    sprich = Array{Int64}(undef,0);
    rich = Array{Int64}(undef,0);
    clock = Array{Float64}(undef,0);
    events = Array{Float64}(undef,0);
    CID = falses(N,maxits);

    freqe = Array{Float64}(undef,maxits);
    freqn = Array{Float64}(undef,maxits);

    mutstep = zeros(Float64,maxits)
    #NOTE strength matrix can't be built a-priori!
    # #Build the strength matrix apriori
    # strength = vec(pi*sum(nb0,dims=2)) .- vec(sqrt(2)*sum(eb,dims=2)) .- vec(sum(eb,dims=1));
    # smatrix = Array{Float64}(copy(eb));
    # for i=1:size(eb)[1]
    #     smatrix[i,:] .= eb[i,:] * strength[i];
    # end
    # # smatrix[findall(iszero,eb)] = NaN;
    #
    minstrength = -ce*Float64(S) - cp*Float64(S);

    evolutiontable = [[0 0 0 1 1 1 2 2 2 3 3 3];[1 2 3 0 2 3 0 1 3 0 1 2]];
    tallytable = [4.1 4.2 5.1 4.3 4.4 5.2 4.5 4.6 5.3 6.1 6.2 6.3];


    t=0;
    it = 0;
    while it < maxits

        #does sorting CID make a difference?
        sort!(cid);

        cid_old = copy(cid);
        #Which are species?
        spcid = intersect(sp_v,cid);
        #Which are objects?
        ocid = setdiff(cid,spcid);


        #COUNT POTENTIAL COLONIZERS
        trophiclinked = setdiff(int_id[(sum(eb[:,[1;cid]],dims=2) .> 0)[:,1]],cid);
        #For each trophiclinked, count number of assimilate and need interactions in system
        #Determine in the proportion that already exists is >= the threshold
        a_fill = ((sum(eb[trophiclinked,[1;cid]],dims=2)./sum(eb[trophiclinked,:],dims=2)) .>= e_t)[:,1];
        prop_n = sum(nb0[trophiclinked,[1;cid]],dims=2)./sum(nb0[trophiclinked,:],dims=2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] .= 1;
        n_fill = (prop_n .>= n_t)[:,1];
        #Build a list of those that pass each individual test
        a_pass = trophiclinked[a_fill];
        n_pass = trophiclinked[n_fill];
        #Build a list of those that pass both tests (potential colonizers)
        col = intersect(a_pass,n_pass);
        #Count the number that pass
        lcol = length(col);


        #COUNT POTENTIAL EXTINCT SPECIES
        #1) By not fulfilling Eat/Need thresholds
        #Re-calculate e_t and n_t (> e_t; >= n_t)
        a_fill = ((sum(eb[spcid,[1;cid]],dims=2)./sum(eb[spcid,:],dims=2)) .> e_t)[:,1];
        prop_n = sum(nb0[spcid,[1;cid]],dims=2)./sum(nb0[spcid,:],dims=2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] .= 1;
        n_fill = (prop_n .>= n_t)[:,1];
        #Build a list of those that pass each individual test
        a_pass = spcid[a_fill];
        n_pass = spcid[n_fill];
        #which species pass?
        survivors = intersect(a_pass,n_pass);
        spext1 = setdiff(spcid,survivors);

        if length(spcid) > 0

            # Build the strength matrix at each community state
            # Strength values change over time so need to ve updated
            # Only record strength values of species (hence the [1:length(spcid)])
            #NOTE: Needs won't change; Eats is based on POTENTIAL niche; Vuln changes per timestep
            strength = vec(cn*sum(nb0[spcid,cid],dims=2)) .- vec(ce*sum(eb[spcid,:],dims=2)) .- (vec(cp*sum(eb[spcid,cid],dims=1))[1:length(spcid)]);

            cmatrix = Array{Float64}(undef,length(spcid),length(cid));
            for i=1:length(spcid)
                cmatrix[i,:] .= eb[spcid[i],cid] * strength[i];
            end

            # cmatrix = eb[spcid,cid] .* reshape(repeat(strength,outer=length(cid)),length(spcid),length(cid));

            # cmatrix = (eb[spcid,cid]' * (Matrix{Float64}(I,length(strength),length(strength)) .* strength))'

            #'zero' entrees need to be lower than any possible strength
            #So they will be effectively ignored
            #-sqrt(2)*S-S is the theoretical min strength
            cmatrix[cmatrix.==0] .= minstrength;

            #Maximum strength values for each resource utilized
            cmax = findmax(cmatrix,dims=1)[1];

            prext_comp = trues(length(spcid));
            for i=1:length(spcid)

                #Don't count sun
                #catalogue prey for all species/objects in the system (excluding sun)
                ieats = Array{Bool}(eb[spcid[i],cid]);
                #If you have >= the max strength for any of those prey, you stay
                #This means that a pure primary producer is not evaluated
                prext_comp[i] = any(ieats)*(any(strength[i] .>= cmax[ieats])==false);

            end
            spext2 = spcid[prext_comp];
        else
            spext2 = Array{Int64}(undef,0);
        end

        primext = copy(spext2);
        secext = copy(spext1);
        #
        # spext2 = Array{Int64}(0);
        spext = unique([spext1;spext2]);
        lspext = length(spext);



        #COUNT POTENTIAL EXTINCT OBJECTS
        obext = ocid[findall(iszero,vec(sum(m_b[spcid,ocid],dims=1)))];
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
        push!(clock,t);

        #Choose a random event
        re = rand();
        tally = -10;

        if re < (lcol/levents)

            #COLONIZATION FUNCTION
            c_sp = rand(col);
            #select made objects as well
            c_ob_full = findall(isodd,m_b[c_sp,:]);
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

        if re > ((lcol + lspext)/levents) && re < ((lcol + lspext + lobext)/levents)

            #OBJECT EXTINCTION FUNCTION
            ob_bye = rand(obext,1);
            cid = setdiff(cid_old,ob_bye);
            #tally
            tally = 3;
        end
        
        if re > ((lcol + lspext + lobext)/levents)

            #SPECIES MUTATION
            #Choose species
            #select ANY species
            # spmut = rand(collect(2:S));
            #select FROM COMMUNITY
            spmut = rand(spcid);
            #choose interaction to mutate
            #Mutate ANY interaction
            # intmut = rand(setdiff(collect(2:S),spmut));
            #OR Mutate realized interactions
            intmut = rand(setdiff(cid,spmut));
            oldint = intm[spmut,intmut];
            #choose new interaction
            if in(intmut,spcid)
                #For species, randomly choose ignore, eat, need,
                newint = rand(setdiff([0,1,2],oldint))[1];
            else
                #For objects, randomly choose ignore, eat, need, make
                newint = rand(setdiff([0,1,2,3],oldint))[1];
            end
            
            intm_mut = copy(intm);
            #Update interaction matrix
            intm_mut[spmut,intmut] = newint;
            ebmut,nbmut,nb0mut = intbool(intm_mut);
            tally = tallytable[(evolutiontable[1,:] .== oldint) .& (evolutiontable[2,:] .== newint)][1];

            
            #Does the mutant outcompete the parent?
            # if in(spmut,spcid)

            spmutloc = findall(x->x==spmut,spcid)[1];
            #Calculate strength
            strength_mut = (cn*sum(nb0mut[spmut,cid])) .- (ce*sum(ebmut[spmut,:])) .- (cp*sum(ebmut[cid,spmut]));


            if strength_mut > strength[spmutloc]
                #accept mutation
                intm = copy(intm_mut);
                eb = copy(ebmut);
                nb = copy(nbmut);
                nb0 = copy(nb0);
                i_b = copy(i_bmut);

                mutstep[it] = 1.0;

            end
                


        end

        freqe[it] = sum(eb[spcid,cid])/(length(cid));
        freqn[it] = sum(nb0[spcid,cid])/(length(cid));

        #NOTE - updating CID....
        CID[cid,it] .= true;
        push!(sprich,length(spcid));
        push!(rich,length(cid));
        push!(events,tally);
    end #end time steps

    intm_evo = copy(intm);



    return(
    sprich,
    rich,
    clock,
    CID,
    intm_evo,
    mutstep,
    freqe,
    freqn,
    events
    )

end
