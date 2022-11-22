"""
    setUpPool(S, lambda, numBasalResources, SSProbs, SMProbs, diverse)

    Returns a random pool network with 'S' species, 'numBasalResources' basal resources and an average of 'lambda' makes per species.
    The probabilities of species-species interactions are given by 'SSProbs' and the probabilities of species modifier and species-basal resource interactions bei 'SMProbs'.
"""
function setUpPool(S, lambda, numBasalResources, SSProbs, SMProbs, diverse)
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
        addSpec!( poolnet, specId, ENIgMaVert() );
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

    # compute number of link types between node types
    numModEats = Int(round(SMProbs.p_e*S*numMods))  #species-modifier eats
    numModNeeds = Int(round(SMProbs.p_n*S*numMods)) #species-modifier needs etc.
    numBasalResEats = Int(round(SMProbs.p_e*S*numBasalResources)) #no needs for basal Resources
    numSpecEats = Int(round(SSProbs.p_e*S*(S-1)))
    numSpecNeeds = Int(round(SSProbs.p_n*S*(S-1)))

    #first create all non make links between species and modifiers
    possibleEdges = [(specId,modId) for specId in poolnet.spec for modId in poolnet.mods]   #set of possible eges
    realizedEdges = sample(possibleEdges, numModEats + numModNeeds, replace = false)    #draw all realizes edges as one sample to avoid double links
    
    for (predatorId, preyId) in view(realizedEdges,1:numModEats)        # make as many species-mod eats as requested
        adde!(poolnet[predatorId],predatorId,poolnet[preyId],preyId)
    end
    for (specId,needId) in view(realizedEdges,(numModEats + 1):(numModNeeds + numModEats))  #make as many species-mod needs as requested
        addn!(poolnet[specId],needId)
    end

    #create all basal resource interactions analogous to above
    possibleEdges = [(specId,basalResId) for specId in poolnet.spec for basalResId in poolnet.basalRes]
    realizedEdges = sample(possibleEdges, numBasalResEats, replace = false)
    for (predatorId, preyId) in view(realizedEdges,1:numModBasalResEats)
        adde!(poolnet[predatorId],predatorId,poolnet[preyId],preyId)
    end

    #create all species-species interactions analogous to above
    possibleEdges = [id1,id2 for id1 in poolnet.spec for id2 in poolnet.spec if id1 != id2]
    realizedEdges = sample(possibleEdges, numSpecEats + numSpecNeeds, replace = false)

    for (predatorId, preyId) in view(realizedEdges,1:numSpecEats)
        adde!(poolnet[predatorId],predatorId,poolnet[preyId],preyId)
    end
    for (specId,needId) in view(realizedEdges,(numSpecEats + 1):(numSpecNeeds + numSpecEats))
        addn!(poolnet[specId],needId)
    end

    #check if all species hav an eat otherwise assign a basal resource as food
    for specId in poolnet.spec
        spec = poolnet[specId]
        if isempty(spec.eat)
            basalEatId = rand(poolnet.basalRes)
            adde!(spec,specId,poolnet[basalEatId],basalEatId)
        end
    end
    return poolnet;
end
