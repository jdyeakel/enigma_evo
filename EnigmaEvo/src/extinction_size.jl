function dictToDistribution(dict)
    minVal,maxVal = extrema(keys(dict))
    offset = -minVal + 1;
    shiftedDistVec = zeros(maxVal + offset)
    numEvents = sum(values(dict))
    for (key,val) in dict
        shiftedDistVec[key + offset] = val/numEvents;
    end
    return OffsetArray(shiftedDistVec,minVal:maxVal)
end

"""
    get_extinction_size_distrib(sprich;break_at_0 = false)

    Computes the distribution of sizes of growth and extinction cascades within the species richness 'sprich'.

    A growth/extinction cascade is considered to be a continuous chain of only extinction/growth events. 
    If 'break_at_zero' = true an event that doesn't affect the species richness ends a cascade. Otherwise they are ignored.
"""
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

    return dictToDistribution(ext_length_dist)
end

"""
    speciesRichnessDeltaPrePostEvents(sprich, itterations, clock, Δt, offset=0)
"""
function deltaSPrePostEventsDist(sprich, itterations, clock, Δt; offset=0)
    ΔSPreDist = Dict{Int,Int}()
    ΔSPostDist = Dict{Int,Int}()
    tMin = Δt + clock[1]
    tMax = clock[end] - Δt
    for (it,t) in ((it,clock[it]) for it in itterations if it > offset )
        if t >= tMin
            start = findfirst(s -> s > t - Δt,clock) - 1
            ΔS = sprich[it] - sprich[start]
            if haskey(ΔSPreDist, ΔS)
                ΔSPreDist[ΔS] += 1
            else
                ΔSPreDist[ΔS] = 1
            end
        end

        if t <= tMax
            windowEnd = findfirst(s -> s > t + Δt, clock) - 1
            ΔS = sprich[windowEnd] - sprich[it]
            if haskey(ΔSPostDist, ΔS)
                ΔSPostDist[ΔS] += 1
            else
                ΔSPostDist[ΔS] = 1
            end
        end
    end

    return dictToDistribution(ΔSPreDist), dictToDistribution(ΔSPostDist)
end


ΔSPrePostEventsDist = deltaSPrePostEventsDist


"""
    speciesRichnessDeltaPrePostEvents(sprich, itterations, clock, Δt, offset=0)
"""
function deltaSPrePostEvents(sprich, itterations, clock, Δt; offset=0)
    ΔSPre = Int[]
    sizehint!(ΔSPre,length(itterations))
    ΔSPost = Int[]
    sizehint!(ΔSPost,length(itterations))
    tMin = Δt + clock[1]
    tMax = clock[end] - Δt
    for (it,t) in ((it,clock[it]) for it in itterations if it > offset )
        if t >= tMin
            start = findfirst(>(t - Δt),clock) - 1
            push!(ΔSPre, sprich[it] - sprich[start])
        end

        if t <= tMax
            windowEndShifted = findfirst(>(t + Δt), clock)
            if windowEndShifted === nothing
                if clock[end] == t + Δt
                    windowEnd = length(clock)
                else
                    continue
                end
            else
                windowEnd = windowEndShifted - 1
            end
            push!(ΔSPost, sprich[windowEnd] - sprich[it])
        end
    end

    return ΔSPre, ΔSPost
end