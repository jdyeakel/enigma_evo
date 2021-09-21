function assemblyevo(edgelist,sID,oID,lambda,
    athresh,nthresh,maxits,probmut,cn,ce,cp)

    S = length(sID);
    N = S + length(oID);
    
    evID = Array{Int64}(undef,0);
    evIDstart = copy(N) + 1; #keep track of evolved species ids - start after objects

    MaxN = convert(Int64,floor(S + S*lambda));
    
    cid = Array{Int64}(undef,0);
    #Add the primary resource
    push!(cid,1);

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
        splist = [sID;evID];
        spcid = intersect(splist,cid);
        #Which are objects?
        ocid = setdiff(cid,spcid);

        # COLONIZATION
        col,lcol = potcol2(edgelist,splist,cid,e_t,n_t);
        # #Species with needs
        # ntest = sort(unique(edgelist[(edgelist[:,1] .<= 200) .& (edgelist[:,3] .== 2),1]))
        # #To see if species with needs are in colonization list (test when cid = 1)
        # indexin(col,ntest)

        # EXTINCTION
        #COUNT POTENTIAL EXTINCT SPECIES ~ PRIMARY
        primext = potextinct2(edgelist,spcid,cid);
        
        #COUNT POTENTIAL EXTINCT SPECIES ~ SECONDARY
        secext = potsecextinct2(edgelist,spcid,cid,cn,ce,cp)

        spext = unique([primext;secext]);
        lspext = length(spext);

        # OBJECTS
        #COUNT POTENTIAL EXTINCT OBJECTS
        obext = ocid[findall(iszero,vec(sum(m_b[spcid,ocid],dims=1)))];
        lobext = length(obext);

        levents = sum([lcol;lspext;lobext]);
        #mutation probability is set (and not dependent on state)
        lmutations = Int64(floor(probmut*levents));
        #No mutations if there are no interactions in the system
        if length(spcid) < 2
            lmutations = 0;
        end
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
            
            # NOTE: Here allow DISRUPTIVE EVOLUTION

            spmutloc = findall(x->x==spmut,spcid)[1];


            # Add new species to the Evolve Matrix - update binary versions

            levolvem = size(evolvem)[1];
            evolvem_new = cat(evolvem,intm_mut[spmut,:]);
            # New species id
            eid += 1;
            espcid = [espcid eid];

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
