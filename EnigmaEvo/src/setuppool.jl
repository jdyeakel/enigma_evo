function setuppool(S, lambda, numBasalResources, SMProbs, SOprobs,diverse) #, OOprobs)
    numBasalResources < 1 && throw(DomainError(numBasalResources, "There has to be at least one basal resource."))
    lambda < 0 && throw(DomainError(lambda, "lambda must be non-negative."))

    numMakes = Int(round(lambda*S));    #num of original makes, also used as the number of possible objects
    maxN = S + numBasalResources + numMakes      # max num of initial vertices
    estsize = maxN + Int(round(maxN*1*diverse)); # give a guess of the max system size might be optimized by tuning this value with observations
    poolnet::ENIgMaGraph = ENIgMaGraph(estsize); # create empty pool network
    
    #create basal resources
    for i in 1:numBasalResources
        basResId = getnextid!(poolnet)
        addBasalRes!( poolnet, basResId, ENIgMaVert());
    end

    #create species, fill later
    #firstSpecId = numBasalResources + 1
    for i in 1:S
        specId = getnextid!(poolnet)
        addspec!( poolnet, specId, ENIgMaVert() );
    end

    potMods = [ENIgMaVert() for i in 1:numMakes];   #create vector of potential modifiers. Add them to poolnet later if they are made by a species

    for i in 1:numMakes                     #now create all the makes
        engineerId = rand(poolnet.spec)
        createdMod = rand(potMods)
        addn!(createdMod,engineerId)        #for now only save the reciprocal need of the modifier as the id of the modifier isnt yet specified
    end

    for potMod in potMods                   # go through all modifiers and....
        if !isEmpty(potMod.need)            # ...check if they are created by any species
            modId = getnextid!(poolnet)     # give the modifier an id
            for engineerId in potMod.need   # add all makes of that modifiers engineers
                addm!(poolnet[engineerId],modId)
            end
            addMod!(poolnet, modId, potMod) # add modifier to pool net
        end
    end

    numMods = getNumMods(poolnet) #number of modifiers

    # compute number of link types between species and (mods + basal res) and spec and spec
    numModBasalResEats = Int(round(SMProbs.p_e*S*(numBasalResources+numMods)))
    numModBasalResNeeds = Int(round(SMProbs.p_n*S*(numBasalResources+numMods)))
    numSpecEats = Int(round(SSProbs.p_e*S*(S-1)))
    numSpecNeeds = Int(round(SSProbs.p_n*S*(S-1)))

    modsAndBasalRes = union(poolnet.mods,poolnet.basalRes)  #set of all modifier and basal resource ids
    for i in 1:numModBasalResEats
        predatorId = rand(poolnet.spec)
        preyId = rand(modsAndBasalRes)
        adde!(poolnet[predatorId],predatorId,poolnet[preyId],preyId)
    end

    for i in 1:numModBasalResNeeds
        specId = rand(poolnet.spec)
        needId = rand(modsAndBasalRes)
        addn!(poolnet[specId],needId)
    end

    for i in 1:numSpecEats
        predatorId,preyId = sample(poolnet.spec,2, replace=false)
        
        adde!(poolnet[predatorId],predatorId,poolnet[preyId],preyId)
    end




    

    funcs = [adde!,addn!,adde!]


    lambdaSS = S*(SSprobs.p_e + SSprobs.p_n);               # Expectet number of (interesting) SO interactions (only in) per species

    pdistSS = Poisson(lambdaSS);
    SSintpS = rand(pdistSS,S);            # draw number of SS interactions for each species

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
            deln!(newspec,1)    #cant have double links
            adde!(newspec,1)    #make autotroph!
            addf!(poolnet[1], specid)
        end
    end

    #prepare first species (add basal resource as eat discard all need interactions to assure that at least one species can colonize)
    autotroph = poolnet[2];
    adde!(autotroph,1);
    addf!(poolnet[1], 2);
    empty!(autotroph.need);


    return poolnet;
end
