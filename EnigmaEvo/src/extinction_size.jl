function get_extinction_size_distrib(sprich,break_at_0 = false)    #pretty shure there is a faster solution...and a shorter but more confusing one and probably slower(using diff and diff again) but
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
        elseif break_at_0
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