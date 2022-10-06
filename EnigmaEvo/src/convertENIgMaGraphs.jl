using .ENIgMaGraphs

function converttointeractionmat(g::ENIgMaGraph)
    N = length(g);
    intm = zeros(Int,N,N);

    for (id,v) in g
        if g.hasspec[id]
            intm[id,id] = 2;
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

    intm[1,1] = 0;

    return intm;
end

function converttoENIgMaGraph(intm)
    N = size(intm)[1];
    g = ENIgMaGraph(N,IdManager(N));

    for id in 1:N
        if intm[id,id] == 2
            addspec!(g,id,ENIgMaVert())
        else
            addmod!(g,id,ENIgMaVert())
        end
    end


    for row in 1:N
        isspec = false;
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

    addn!(g[1],1);
    return g;
end