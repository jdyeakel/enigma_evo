function setuppool(S, lambda, SSprobs, SOprobs) #, OOprobs)
    # first draw all SO connections, especially all make connections, all objects not created by any species will be discarded
    numOpos = Int(round(lambda*S));    #num obj possible

    lambdaSS = S*(SSprobs.p_e + SSprobs.p_n);               # Expectet number of (interesting) SO interactions (only in) per species
    p_m = lambda/numOpos;  
    
    probSOint = SOprobs.p_e + SOprobs.p_n + p_m     #probability of any SO interaction (can be added because they are mutually exclusive)

    lambdaSO = numOpos*probSOint;      # Expectet number of (interesting) SO interactions (only in) per species
    intprobs = [SOprobs.p_e,SOprobs.p_n,p_m]/probSOint;     #conditional (link exists) probabilities of interaction types
    intprobthresh = cumsum(intprobs);   #used as thresholds

    #A species is an engineer if number of objects > 0
    pdistSS = Poisson(lambdaSS);
    pdistSO = Poisson(lambdaSO);
    
    SSintpS = rand(pdistSS,S);            # draw number of SO interactions for each species
    SOintpS = rand(pdistSO,S);            # draw number of SS interactions for each species

    objectrealized = falses(numOpos);           # checks wheather any species creates an object
    SOintmat = zeros(Int,numOpos,S);            # faster if fast index (innermost) goes along coloumns (changes rows) in julia...
    #obline = collect(Int64,1:numOpos); actually slower than using 1:numOpos directly
    for (specid,numOSint) = enumerate(SOintpS)          #loop over all species and their numbers of SO interactions
        if numOSint > 0                                 #has any interactions?
            objids = sample(1:numOpos,numOSint,replace=false);   #draw objects to interact with
            for objid in objids                         #loop over these
                inttype = rand()                        #draw number that decides the type of interaction
                if inttype < intprobthresh[1]           #eat interaction?
                    SOintmat[objid,specid] = 1;
                elseif inttype < intprobthresh[2]       #need interaction?
                    SOintmat[objid,specid] = 2;         
                else                                    #make interaction
                    SOintmat[objid,specid] = 3;
                    objectrealized[objid] = true;       #note that object is made by a species
                end
            end
        end
    end

    # discard species that are not made by any species (note that this does not change the average number of made interactions per species
    # and neither the expected probability of need and eat SO links)
    SOintmat = SOintmat[objectrealized,:]

    O = size(SOintmat,1);   #num objects
    
    N = 1 + S + O;  # basal resource neither Species nor object here
    estsize = N + Int(round(N*.1*diverse)); #give a guess of the max system size might be optimized by tuning this value with observations
    poolnet::ENIgMaGraph = ENIgMaGraph(estsize);    # create empty pool network

    #create basal resource
    basalres = ENIgMaVert()
    #add itself as need in order to  make sure it doesnt go extinct (needs one need as object...easier to program that way)
    addn!(basalres,1);
    addmod!( poolnet, 1, basalres);
    


    #create species, fill later
    for specid in 2:(S+1)
        addspec!( poolnet, specid, ENIgMaVert() );
    end

    #Create modifieres/objects, fill later
    for objid in (S + 2):N
        addmod!( poolnet, objid, ENIgMaVert() );
    end

    for (shiftedspecid,numSSint) = enumerate(SSintpS)   #loop over all species and their numbers of SS interactions
        specid = shiftedspecid + 1;
        newspec = poolnet[specid]
        if numSSint > 0                                 #has any interactions?
            intids = sample(1:S,numSSint,replace=false);         #draw species or basal resource to interact with (should be each time over different set to avoid self interactions, taken care of later)
            for intid in intids                         #loop over these
                if intid >= specid                      #convert index to avoid self interactions (note that proper spec ids are from 2:(S+1))
                    intid +=1;
                end
                inttype = rand()                        #draw number that decides the type of interaction
                if inttype < SSprobs.p_e/(SSprobs.p_e + SSprobs.p_n)    #eat interaction?
                    adde!(newspec,intid)
                    addf!(poolnet[intid], specid);
                else                                     #need interaction?
                    addn!(newspec,intid)
                end
            end
        end
        if SOintpS[shiftedspecid] > 0   #might be optimized by deleting this if
            for (shiftedobjid,inttype) in enumerate(SOintmat[:,shiftedspecid])
                objid = shiftedobjid + S + 1;
                if inttype == 1
                    adde!(newspec,objid)
                    addf!(poolnet[objid],specid)
                elseif inttype == 2
                    addn!(newspec, objid)
                elseif inttype == 3
                    addm!(newspec,objid)
                    addn!(poolnet[objid],specid)
                end
            end
        end
        if isempty(newspec.eat) #no eat interaction?
            adde!(newspec,1)    #make autotroph!
            addf!(poolnet[1], specid)
        end
    end

    #prepare first species (add basal resource as eat discard all need interactions to assure that at least one species can colonize)
    autotroph = poolnet[2];
    adde!(autotroph,1);
    addf!(poolnet[1], 2);
    empty!(autotroph.need);

    ENIgMaGraphs.setmaxid!(poolnet.idmanager,N);    # dont like the design so far

    return poolnet;
end
