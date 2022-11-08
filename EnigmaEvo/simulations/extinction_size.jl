function get_extinction_size_distrib(sprich)    #pretty shure there is a faster solution...and a shorter but more confusing one and probably slower(using diff and diff again) but
    del_sprich = diff(sprich);
    ext_length_dist = Dict{Int,Int}()
    increase = 0;
    decrease = 0;
    for delta in del_sprich
        if delta > 0
            if decrease < 0
                if haskey(ext_length_dist,decrease)
                    ext_length_dist[decrease] += 1
                else
                    ext_length_dist[decrease] = 1
                end
                    decrease = 0
            end
            increase += delta
        elseif delta < 0
            if increase > 0
                if haskey(ext_length_dist,increase)
                    ext_length_dist[increase] += 1
                else
                    ext_length_dist[increase] = 1
                end
                increase = 0
            end
            decrease += delta
        else
            if decrease < 0
                if haskey(ext_length_dist,decrease)
                    ext_length_dist[decrease] += 1
                else
                    ext_length_dist[decrease] = 1
                end
                decrease = 0
            end
            if increase > 0
                if haskey(ext_length_dist,increase)
                    ext_length_dist[increase] += 1
                else
                    ext_length_dist[increase] = 1
                end
                increase = 0
            end
            if haskey(ext_length_dist,0)
                ext_length_dist[0] += 1
            else
                ext_length_dist[0] = 1
            end
        end
    end
    if decrease < 0
        if haskey(ext_length_dist,decrease)
            ext_length_dist[decrease] += 1
        else
            ext_length_dist[decrease] = 1
        end
        decrease = 0
    end
    if increase > 0
        if haskey(ext_length_dist,increase)
            ext_length_dist[increase] += 1
        else
            ext_length_dist[increase] = 1
        end
        increase = 0
    end
    max_decr,max_incr = extrema(keys(ext_length_dist))
    offset = -max_decr + 1;
    ext_len_dist_vec = zeros(max_incr - max_decr + 1)
    num_events = sum(values(ext_length_dist))
    for (key,val) in ext_length_dist
        ext_len_dist_vec[key + offset] = val/num_events;
    end
    return OffsetArray(ext_len_dist_vec,max_decr:max_incr)
end

function get_extinction_size_distrib_test()
    working = true;
    function isDist(dist,caseNum)
        if sum(dist) != 1
            println("Error: extinction_size_distrib does not return probability distribution in test case $(caseNum).")
            println("\tsum(dist) != 1")
            return false
        end
        return true
    end
    dist = extinction_size_distrib(zeros(10))
    working = isDist(dist,1)
    if dist[0] != 1
        working = false;
        println("Error: extinction_size_distrib failed test case 1.")
    end
    dist = extinction_size_distrib(1:10)
    working = isDist(dist,2)
    if dist[9] != 1
        working = false;
        println("Error: extinction_size_distrib failed test case 2.")
    end
    dist = extinction_size_distrib([7,7,6,5,4,3,3,4,5,4])
    working = isDist(dist,3)
    if dist[[-4,-1,0,2]] != [1,1,2,1.]/5
        working = false;
        println("Error: extinction_size_distrib failed test case 3.")
    end
    return working
    dist = extinction_size_distrib([1,2,3,2,3,4,3,2,3,4,5,6,5,4,5,4,3])
    working = isDist(dist,3)
    if dist[[-2,-1,1,2,4]] != [3,1,1,2,1]/8
        working = false;
        println("Error: extinction_size_distrib failed test case 3.")
    end
    return working
end