function assemblyevo_diverse(edgelist,sID,oID,e_t,n_t,maxits,probmut,cn,ce,cp,div_t)

    S = length(sID);
    N = S + length(oID);
    
    evID = Array{Int64}(undef,0);
    evIDticker = copy(N) + 2; #keep track of evolved species ids - start after objects

    # MaxN = convert(Int64,floor(S + S*lambda));
    
    cid = Array{Int64}(undef,0);
    #Add the primary resource
    push!(cid,1);

    # clock = Array{Float64}(undef,0);
    # CID = falses(N,maxits);

    sprich = Array{Int64}(undef,maxits);
    obrich = Array{Int64}(undef,maxits);
    events = Array{Float64}(undef,maxits);
    clock = Array{Float64}(undef,maxits);
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
    # minstrength = -ce*Float64(S) - cp*Float64(S);
    evolutiontable = [[0 0 0 1 1 1 2 2 2 3 3 3];[1 2 3 0 2 3 0 1 3 0 1 2]];
    tallytable = [4.1 4.2 5.1 4.3 4.4 5.2 4.5 4.6 5.3 6.1 6.2 6.3];

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
        ocid = setdiff(cid,[1;spcid]);

        # COLONIZATION
        col = potcol2(edgelist,splist,cid,e_t,n_t);
        lcol = length(col);
        # #Species with needs
        # ntest = sort(unique(edgelist[(edgelist[:,1] .<= 200) .& (edgelist[:,3] .== 2),1]))
        # #To see if species with needs are in colonization list (test when cid = 1)
        # indexin(col,ntest)

        # EXTINCTION
        #COUNT POTENTIAL EXTINCT SPECIES ~ PRIMARY
        primext = potextinct2(edgelist,spcid,cid);
        
        #COUNT POTENTIAL EXTINCT SPECIES ~ SECONDARY
        secext = potsecextinct2(edgelist,spcid,cid,cn,ce,cp);

        spext = unique([primext;secext]);
        lspext = length(spext);

        # OBJECTS
        #COUNT POTENTIAL EXTINCT OBJECTS
        obext = potobextinct2(edgelist,ocid,cid);
        lobext = length(obext);

        #Possible ecological events: colonization, sp exitinction, ob extinction
        lecoevents = sum([lcol;lspext;lobext]);

        #No mutations if there are no interactions in the system
        if length(spcid) < 2
            probevo = 0.;
        else
            probevo = copy(probmut);
        end
        
        #Total number of events include ECO events + EVO event
        levents = Int64(floor((lecoevents)/(1-probevo))); 
        # lmutations = levents - lecoevents;
        
        
        
        #Redefine levents to include lmutations
        # levents = sum([lcol;lspext;lobext;lmutations]);

        

        #Choose a random event
        re = rand();
        tally = -10;

        if re < (lcol/levents)

            #COLONIZATION FUNCTION
            c_sp = rand(col);
            #select made objects as well
            c_ob_full = intfind(edgelist,c_sp,3);
            # c_ob_full = findall(isodd,m_b[c_sp,:]);
            #Only bring in objects that aren't already in the community
            c_ob = setdiff(c_ob_full,cid_old);
            colonize = [c_sp;c_ob];
            #Update community
            cid = [cid_old;colonize];
            #Update species (we do this above)
            # spcid = [spcid;c_sp];
            #tally
            tally = 0.;

            # dt = 1/levents;
            # t += dt;
            # it += 1;
            # push!(clock,t);
        end

        if (re > (lcol/levents)) && (re < ((lcol + lspext)/levents))

            #SPECIES EXTINCTION FUNCTION
            #select species to go extinct
            sp_bye = rand(spext,1);
            cid = setdiff(cid_old,sp_bye);
            # spcid = setdiff(spcid,sp_bye);
            #tally
            if in(sp_bye[1],primext)
                #These are extinctions from competitive exclusion
                tally = 1.;
            end
            if in(sp_bye[1],secext)
                #These are secondary extinctions
                tally = 2.;
            end
            # dt = 1/levents;
            # t += dt;
            # it += 1;
            # push!(clock,t);
        end

        if (re > ((lcol + lspext)/levents)) && (re < ((lcol + lspext + lobext)/levents))

            #OBJECT EXTINCTION FUNCTION
            ob_bye = rand(obext,1);
            cid = setdiff(cid_old,ob_bye);
            #tally
            tally = 3.;
            # dt = 1/levents;
            # t += dt;
            # it += 1;
            # push!(clock,t);
        end
        
        if (re > ((lcol + lspext + lobext)/levents))

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
            #NOTE: we can't choose from edgelist, because we want to allow non-ineractions to evolve as well
            intmut = rand(setdiff(cid,[1,spmut]),1)[1];
            oldintpos, oldint = intfind_inout(edgelist,spmut,intmut);
            if in(intmut,spcid)
                #For species, randomly choose ignore, eat, need,
                newint = rand(setdiff([0,1,2],oldint))[1];
            else
                #For objects, randomly choose ignore, eat, need, make
                newint = rand(setdiff([0,1,2,3],oldint))[1];
            end


            #NOTE: Mutate the primary linneage to a daughter linneage
            spint = edgelist[(edgelist[:,1].==spmut),:];
            spint_rev = edgelist[(edgelist[:,2].==spmut),:];
            spmut_daughter = copy(evIDticker);
            #Set interactions for new daughter species
            spint[:,1] .= spmut_daughter;
            spint_rev[:,2] .= spmut_daughter;
            evIDticker += 1;
            newedge = [spmut_daughter intmut newint];
            #If interaction did not exist, add it
            if oldint == 0
                spint_daughter = vcat(spint,newedge);
            else
            #If interaction does exist, mutate it
                mutintpos, mutint = intfind_inout(spint,spmut_daughter,intmut)
                spint_daughter = copy(spint);
                spint_daughter[mutintpos,3] = newint;
                #Delete non-interaction if it exists
                spint_daughter = spint_daughter[(spint_daughter[:,3] .!= 0),:];
                # If we have mutated a make interaction, we need to 'turn off' the reverse need interaction
                if oldint == 3
                    mutintpos_rev, mutint_rev = intfind_inout(spint_rev,intmut,spmut_daughter)
                    spint_rev = spint_rev[setdiff(1:size(spint_rev)[1],mutintpos_rev),:]
                end
            end

            daughter_interactions = [spint_daughter;spint_rev];
            
            #Update edgelist
            edgelist = vcat(edgelist,daughter_interactions);
            #Update cid
            push!(cid,spmut_daughter);
            #Update evID
            push!(evID,spmut_daughter);
            #Update oID?
            #NOTE: It would be good to 'evolve' objects, but we don't have this feature yet

            #Record the interaction evolutionary transition
            tally = tallytable[(evolutiontable[1,:] .== oldint) .& (evolutiontable[2,:] .==newint)][1];

            #Does the parent linneage survive? NOTE: This is permanent extinction
            parent_strength = strength(edgelist,spmut,cid,cn,ce,cp);
            daughter_strength = strength(edgelist,spmut_daughter,cid,cn,ce,cp);
            rdiversify = rand();
            if rdiversify > div_t
                #Eliminate lesser competitor
                if parent_strength < daughter_strength
                    intpos = findall(isone,(edgelist[:,1] .== spmut) .| (edgelist[:,2] .== spmut));
                    edgelist = edgelist[setdiff(1:size(edgelist)[1],intpos),:];
                    #Remove spmut from cid
                    cid = setdiff(cid,spmut);
                    #Remove from either sID or evID
                    sID = setdiff(sID,spmut);
                    evID = setdiff(evID,spmut);
                else
                    intpos = findall(isone,(edgelist[:,1] .== spmut_daughter) .| (edgelist[:,2] .== spmut_daughter));
                    edgelist = edgelist[setdiff(1:size(edgelist)[1],intpos),:];
                    #Remove spmut from cid
                    cid = setdiff(cid,spmut_daughter);
                    #Remove from either sID or evID
                    sID = setdiff(sID,spmut_daughter);
                    evID = setdiff(evID,spmut_daughter);
                end
            end
            # If rdiversity < div_t, then keep both
            mutstep[it] = 1.0;
            
            
            
        end

        dt = 1/levents;
        t += dt;
        it += 1;
        events[it] = tally;
        clock[it] = t;

        #richness
        sprich[it] = length(intersect(cid,[sID;evID]));
        obrich[it] = length(intersect(cid,oID));
        #Link density
        freqe[it] = sum(edgelist[:,3] .== 1)/length(cid);
        freqn[it] = sum(edgelist[:,3] .== 2)/length(cid);

        #NOTE - updating CID....
        # CID[cid,it] .= true;
       
    end #end time steps


    return(
    sprich,
    obrich,
    clock,
    edgelist,
    cid,
    evID,
    mutstep,
    freqe,
    freqn,
    events
    )

end
