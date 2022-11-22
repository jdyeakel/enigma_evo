using .ENIgMaGraphs

"""
    converttointeractionmat(g::ENIgMaGraph)

    Returns 'g' in matrix form.
"""
function converttointeractionmat(g::ENIgMaGraph)
    N = g.idmanager.maxid;
    intm = zeros(Int,N,N);

    for (id,v) in g
        if g.hasBasalRes[id]
            intm[id,id] = 1;
        elseif g.hasspec[id]
            intm[id,id] = 2;
        elseif g.hasMod[id]
            intm[id,id] = 3
        end

        for eid in v.eat
            intm[id,eid] = 1;
        end
        for nid in v.need
            intm[id,nid] = 2;
        end
        for mid in v.make
            intm[id,mid] = 3;
        end
    end

    return intm;
end

"""
    converttoENIgMaGraph(intm)

    Converts the interaction matrix 'intm' to an ENIgMaGraph and returns it.
"""
function converttoENIgMaGraph(intm,estSize=nothing)
    N = size(intm)[1];
    g = ENIgMaGraph(estSize === nothing ? N : estSize);

    for i in 1:N
        vertType = intm[i,i]
        id = getNextId!(g)  #get next id even if species already extinct to hopefully get inverse of converttointeractionmat
        if vertType == 2
            addSpec!(g,id,ENIgMaVert())
        elseif vertType == 1 
            addBasalRes!(g,id,ENIgMaVert())
        elseif vertType == 3
            addMod!(g,id,ENIgMaVert())
        end
    end


    for row in 1:N
        if g.hasspec[row]
            newv = g[row];
            for col in 1:N
                if row != col
                    int = intm[row,col];
                    if int == 1
                        adde!(newv, col);
                        addf!(g[col],row);
                    elseif int == 2
                        addn!(newv,col);
                    elseif int == 3
                        addm!(newv,col);
                        addn!(g[col],row);
                    end
                end
            end
        end
    end
    return g;
end