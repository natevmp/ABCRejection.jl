function getParticleDistancesPerMetric(
        distanceFunc::Function,
        particle_tid::AbstractVector,
        dataMetrics,
        nMetrics::Int;
    )
    distance_tid_mid = Array{Float64,2}(undef, length(particle_tid), nMetrics)
    #! Potential change: infer distance type from distanceFunc
    @debug "Computing particle distances" nMetrics
    for (tid, particle) in enumerate(particle_tid)
        distance_mid = distanceFunc(particle.simResults, dataMetrics)
        
        nReturned = distance_mid isa Number ? 1 : length(distance_mid)
        if nReturned != nMetrics
            throw(DimensionMismatch(
                "Function $distanceFunc returned $nReturned distances for particle $tid; expected $nMetrics."
            ))
        end
        distance_tid_mid[tid, :] .= distance_mid
    end
    return distance_tid_mid
end

"""
    competitionRanks!(rank_tid, distance_tid, tid_tidSorted)

Write standard competition ranks for one metric into `rank_tid`, reusing
`tid_tidSorted` as the sorting buffer. All vectors must have matching one-based axes.
Distances must not contain `NaN`. Signed zeros share a rank; infinities are allowed.
The distances are not modified. Return `rank_tid`.
"""
function competitionRanks!(
        rank_tid::AbstractVector{<:Integer},
        distance_tid::AbstractVector,
        tid_tidSorted::AbstractVector{<:Integer},
    )
    sortperm!(tid_tidSorted, distance_tid)

    rank = 1
    for position in eachindex(tid_tidSorted)
        tid = tid_tidSorted[position]

        if position > 1
            previousTid = tid_tidSorted[position - 1]
            if distance_tid[previousTid] < distance_tid[tid]
                rank = position
            end
        end

        rank_tid[tid] = rank
    end
    return rank_tid
end

"""
    sortByMaxRank!(rank_tid_mid)

Return particle indices ordered by their ranks from worst to best. Compare
the largest rank first, then the next-largest rank, until a difference is found.
Lower ranks are better. Identical sorted rank profiles retain input order.
Each row of the input matrix is sorted in place from largest to smallest rank.
"""
function sortByMaxRank!(rank_tid_mid::AbstractMatrix{<:Integer})
    # Reuse the matrix; its columns now represent worst-to-best rank positions.
    rank_tid_rid = rank_tid_mid
    sort!(rank_tid_rid; dims=2, rev=true)

    function lessByWorstRanks(elem1, elem2)
        for rid in axes(rank_tid_rid, 2)
            elem1Rank = rank_tid_rid[elem1, rid]
            elem2Rank = rank_tid_rid[elem2, rid]
            if elem1Rank != elem2Rank
                return elem1Rank < elem2Rank
            end
        end
        return false
    end

    @debug "Particle ranks from worst to best" rank_tid_rid
    return sortperm(axes(rank_tid_rid, 1); lt=lessByWorstRanks)
end

"""
    rankByMax(distance_tid_mid)

Return particle indices ordered by worst-to-best competition rank profiles.
Use one rank matrix, sorting its rows in place for tie-breaking. The distance
matrix is not modified.
"""
function rankByMax(distance_tid_mid::AbstractMatrix)
    nParticles, nMetrics = size(distance_tid_mid)
    rank_tid_mid = Matrix{Int}(undef, nParticles, nMetrics)
    tid_tidSorted = Vector{Int}(undef, nParticles)

    for mid in axes(distance_tid_mid, 2)
        competitionRanks!(
            @view(rank_tid_mid[:, mid]),
            @view(distance_tid_mid[:, mid]),
            tid_tidSorted,
        )
    end

    # Snapshot only when debug logging is enabled, before the in-place sort.
    @debug "Particle ranks by metric" rank_tid_mid=copy(rank_tid_mid)
    return sortByMaxRank!(rank_tid_mid)
end

"""
    rankBySum(distance_tid_mid)

Return particle indices ordered by summed competition ranks. Equal sums retain
input order. Accumulate one metric at a time using reusable rank and sorting
buffers; no full rank matrix is created. The distance matrix is not modified.
"""
function rankBySum(distance_tid_mid::AbstractMatrix)
    nParticles = size(distance_tid_mid, 1)
    score_tid = zeros(Int, nParticles)
    rank_tid = Vector{Int}(undef, nParticles)
    tid_tidSorted = Vector{Int}(undef, nParticles)

    for mid in axes(distance_tid_mid, 2)
        competitionRanks!(
            rank_tid,
            @view(distance_tid_mid[:, mid]),
            tid_tidSorted,
        )
        score_tid .+= rank_tid
    end

    @debug "Summed ranks per particle" score_tid
    return sortperm(score_tid)
end

"""
    rankParticles(
        distanceFunc::Function,
        particle_tid::AbstractVector,
        dataMetrics;
        nMetrics::Union{Nothing,Int}=nothing,
        algorithm::Symbol=:max,
    )

Return particle indices ordered by their per-metric ranks, with lower ranks better.
Each metric uses standard competition ranks, with equal distances sharing a rank.
Use `algorithm=:max` (default) to compare each particle's worst rank, breaking
ties by its second-worst rank, then third-worst, and so on. Identical sorted
rank profiles retain input order. Use `algorithm=:sum` to compare the sum of
the ranks; equal sums retain input order without a secondary tie-breaker.

Signed zeros share a rank, and infinite distances are allowed. A particle with
`NaN` in any distance is excluded before ranks are computed. One warning lists
the excluded particle positions. The returned ordering contains only eligible
particles and is empty if none remain.

`particle_tid` may be a vector or a vector view. Returned indices are positions
within the supplied collection, not indices into a view's parent array.

`distanceFunc(particle.simResults, dataMetrics)` measures the distances between
one particle's simulation results and the observed data. `dataMetrics` is passed
unchanged to this function and may be a tuple, vector, or another representation
it understands. If `nMetrics` is omitted, it defaults to `length(dataMetrics)`;
otherwise, `dataMetrics` need not support `length`.

"""
function rankParticles(
        distanceFunc::Function,
        particle_tid::AbstractVector,
        dataMetrics;
        nMetrics::Union{Nothing,Int}=nothing,
        algorithm::Symbol=:max,
    )
    if algorithm !== :max && algorithm !== :sum
        throw(ArgumentError("Unknown ranking algorithm: $algorithm. Use :max or :sum."))
    end

    if isnothing(nMetrics) nMetrics=length(dataMetrics) end

    distance_tid_mid = getParticleDistancesPerMetric(distanceFunc, particle_tid, dataMetrics, nMetrics)

    @debug "Particle distances by metric" distance_tid_mid

    # check for NaNs
    invalidParticle_tid = vec(any(isnan, distance_tid_mid; dims=2))
    tid_valid = nothing
    if any(invalidParticle_tid)
        @warn "Excluding particles with NaN distances" particleIndices=findall(invalidParticle_tid)
        tid_valid = findall(!, invalidParticle_tid)
        isempty(tid_valid) && return tid_valid
        distance_tid_mid = distance_tid_mid[tid_valid, :]
    end

    tid_rankJoint = if nMetrics == 1
        # Both algorithms reduce to distance ordering; `<` keeps signed zeros tied.
        sortperm(@view(distance_tid_mid[:, 1]); lt=<)
    elseif algorithm === :max
        rankByMax(distance_tid_mid)
    elseif algorithm === :sum
        rankBySum(distance_tid_mid)
    end

    if !isnothing(tid_valid)
        tid_rankJoint = tid_valid[tid_rankJoint]
    end

    @debug "Final particle ordering" tid_rankJoint
    return tid_rankJoint
end
