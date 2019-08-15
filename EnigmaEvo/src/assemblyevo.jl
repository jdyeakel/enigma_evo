function assemblyevo(intm,e_b,n_b,i_b,m_b,n_b0,sp_v,int_id,lambda,
    athresh,nthresh,maxits,probmut,cn,ce,cp)

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
    # strength = vec(pi*sum(n_b0,dims=2)) .- vec(sqrt(2)*sum(e_b,dims=2)) .- vec(sum(e_b,dims=1));
    # smatrix = Array{Float64}(copy(e_b));
    # for i=1:size(e_b)[1]
    #     smatrix[i,:] .= e_b[i,:] * strength[i];
    # end
    # # smatrix[findall(iszero,e_b)] = NaN;
    #
    minstrength = -ce*Float64(S) - cp*Float64(S);


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
        trophiclinked = setdiff(int_id[(sum(e_b[:,[1;cid]],dims=2) .> 0)[:,1]],cid);
        #For each trophiclinked, count number of assimilate and need interactions in system
        #Determine in the proportion that already exists is >= the threshold
        a_fill = ((sum(e_b[trophiclinked,[1;cid]],dims=2)./sum(e_b[trophiclinked,:],dims=2)) .>= athresh)[:,1];
        prop_n = sum(n_b0[trophiclinked,[1;cid]],dims=2)./sum(n_b0[trophiclinked,:],dims=2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] .= 1;
        n_fill = (prop_n .>= nthresh)[:,1];
        #Build a list of those that pass each individual test
        a_pass = trophiclinked[a_fill];
        n_pass = trophiclinked[n_fill];
        #Build a list of those that pass both tests (potential colonizers)
        col = intersect(a_pass,n_pass);
        #Count the number that pass
        lcol = length(col);


        #COUNT POTENTIAL EXTINCT SPECIES
        #1) By not fulfilling Eat/Need thresholds
        #Re-calculate athresh and nthresh (> athresh; >= nthresh)
        a_fill = ((sum(e_b[spcid,[1;cid]],dims=2)./sum(e_b[spcid,:],dims=2)) .> athresh)[:,1];
        prop_n = sum(n_b0[spcid,[1;cid]],dims=2)./sum(n_b0[spcid,:],dims=2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] .= 1;
        n_fill = (prop_n .>= nthresh)[:,1];
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
            strength = vec(cn*sum(n_b0[spcid,cid],dims=2)) .- vec(ce*sum(e_b[spcid,:],dims=2)) .- (vec(cp*sum(e_b[spcid,cid],dims=1))[1:length(spcid)]);

            cmatrix = Array{Float64}(undef,length(spcid),length(cid));
            for i=1:length(spcid)
                cmatrix[i,:] .= e_b[spcid[i],cid] * strength[i];
            end

            # cmatrix = e_b[spcid,cid] .* reshape(repeat(strength,outer=length(cid)),length(spcid),length(cid));

            # cmatrix = (e_b[spcid,cid]' * (Matrix{Float64}(I,length(strength),length(strength)) .* strength))'

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
                ieats = Array{Bool}(e_b[spcid[i],cid]);
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

        levents = sum([lcol;lspext;lobext]);
        #mutation probability is set (and not dependent on state)
        lmutations = Int64(floor(probmut*levents));
        #Redefine levents to include lmutations
        levents = sum([lcol;lspext;lobext;lmutations]);

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
            spmut = rand(collect(2:S));
            #choose interaction to mutate
            intmut = rand(setdiff(collect(2:S),spmut));
            oldint = intm[spmut,intmut];
            #choose new interaction
            newint = rand(setdiff(['e','i','n'],oldint)); #intmut

            intm_mut = copy(intm);
            e_bmut = copy(e_b);
            n_bmut = copy(n_b);
            n_b0mut = copy(n_b0);
            i_bmut = copy(i_b);

            # tp_mmut = copy(tp_m);
            # tind_mmut = copy(tind_m);

            #Update interaction matrix
            intm_mut[spmut,intmut] = newint;
            #Update e_b, n_b, i_b
            if oldint == 'i' && newint == 'n'
                i_bmut[spmut,intmut] = 0;
                n_bmut[spmut,intmut] = 1;
                #tally
                tally = 4.1;
            end
            if oldint == 'i' && newint == 'e'
                i_bmut[spmut,intmut] = 0;
                e_bmut[spmut,intmut] = 1;

                # tp_mmut[spmut,intmut] = 1;
                # tind_mmut[spmut,intmut] = 1;
                #tally
                tally = 4.2;
            end
            if oldint == 'n' && newint == 'i'
                n_bmut[spmut,intmut] = 0;
                n_b0mut[spmut,intmut] = 0;
                i_bmut[spmut,intmut] = 1;
                #tally
                tally = 4.3;
            end
            if oldint == 'n' && newint == 'e'
                n_bmut[spmut,intmut] = 0;
                n_b0mut[spmut,intmut] = 0;
                e_bmut[spmut,intmut] = 1;

                # tp_mmut[spmut,intmut] = 1;
                # tind_mmut[spmut,intmut] = 1;
                #tally
                tally = 4.4;
            end
            if oldint == 'e' && newint == 'n'
                e_bmut[spmut,intmut] = 0;
                n_bmut[spmut,intmut] = 1;
                n_b0mut[spmut,intmut] = 1;

                # tp_mmut[spmut,intmut] = 0;
                # tind_mmut[spmut,intmut] = 0;
                #tally
                tally = 4.5;
            end
            if oldint == 'e' && newint == 'i'
                e_bmut[spmut,intmut] = 0;
                i_bmut[spmut,intmut] = 1;

                # tp_mmut[spmut,intmut] = 0;
                # tind_mmut[spmut,intmut] = 0;
                #tally
                tally = 4.6;
            end
            
            #Does the mutant outcompete the parent?
            
            if in(spmut,spcid)

                spmutloc = findall(x->x==spmut,spcid)[1];
                #Calculate strength
                strength_mut = (cn*sum(n_b0mut[spmut,cid])) .- (ce*sum(e_bmut[spmut,:])) .- (cp*sum(e_bmut[cid,spmut]));


                if strength_mut > strength[spmutloc]
                    #accept mutation
                    intm = copy(intm_mut);
                    e_b = copy(e_bmut);
                    n_b = copy(n_bmut);
                    n_b0 = copy(n_b0);
                    i_b = copy(i_bmut);

                    mutstep[it] = 1.0;


                else
                    #With some probability, accept mutation
                    rr = rand()
                    if rr < exp(-4*abs(strength[spmutloc] - strength_mut))
                        #accept mutation
                        intm = copy(intm_mut);
                        e_b = copy(e_bmut);
                        n_b = copy(n_bmut);
                        n_b0 = copy(n_b0);
                        i_b = copy(i_bmut);

                        mutstep[it] = 1.0;
                    end
                end
            else
                #accept mutation
                intm = copy(intm_mut);
                e_b = copy(e_bmut);
                n_b = copy(n_bmut);
                n_b0 = copy(n_b0);
                i_b = copy(i_bmut);

                mutstep[it] = 1.0;
            end




        end

        freqe[it] = sum(e_b[spcid,cid])/(length(cid));
        freqn[it] = sum(n_b0[spcid,cid])/(length(cid));

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
