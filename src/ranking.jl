function getParticleDistancesPerMetric(
        distDataVSim::Function,
        particle_tid::AbstractVector,
        dataMetrics,
        nMetrics::Int;
    )
    distance_tid_mid = Array{Float64,2}(undef, length(particle_tid), nMetrics)
    @debug "Computing particle distances" nMetrics
    for (tid, particle) in enumerate(particle_tid)
        distance_mid = distDataVSim(particle.simResults, dataMetrics)
        
        nReturned = distance_mid isa Number ? 1 : length(distance_mid)
        if nReturned != nMetrics
            throw(DimensionMismatch(
                "Function $distDataVSim returned $nReturned distances for particle $tid; expected $nMetrics."
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
            if isless(distance_tid[previousTid], distance_tid[tid])
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
        distDataVSim::Function,
        particle_tid::Vector,
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

`distDataVSim(particle.simResults, dataMetrics)`: Function to measure the distance between a single particle and the data. It must take two arguments, the first being the the `simResults` saved in each `particle` and the second the being all metrics of the data to which the particle is being compared.

"""
function rankParticles(
        distDataVSim::Function,
        particle_tid::Vector,
        dataMetrics::Tuple;
        nMetrics::Union{Nothing,Int}=nothing,
        algorithm::Symbol=:max,
    )
    if algorithm !== :max && algorithm !== :sum
        throw(ArgumentError("Unknown ranking algorithm: $algorithm. Use :max or :sum."))
    end

    if isnothing(nMetrics) nMetrics=length(dataMetrics) end

    distance_tid_mid = getParticleDistancesPerMetric(distDataVSim, particle_tid, dataMetrics, nMetrics)

    @debug "Particle distances by metric" distance_tid_mid

    tid_rankJoint = if algorithm === :max
        rankByMax(distance_tid_mid)
    elseif algorithm === :sum
        rankBySum(distance_tid_mid)
    end

    @debug "Final particle ordering" tid_rankJoint
    return tid_rankJoint
end

# function testABCParticlesQuantiles(
#         distDataVSim::Function,
#         particle_tid::Vector,
#         dataMetrics,
#         q::Real,
#         constraintsF::Union{Nothing,Function}=nothing;
#     )
#     nMetrics=length(dataMetrics)
#     distance_tid_mid = getParticleDistancesPerMetric(distDataVSim, particle_tid, dataMetrics)
#     q_mid = [quantile((@view distance_tid_mid[:,mid]), q) for mid in 1:nMetrics]
#     accepted_tid = [all(distance_tid_mid[tid,:].<q_mid) for tid in eachindex(particle_tid)]
#     if isnothing(constraintsF)
#         return accepted_tid
#     end
#     acceptedConstraints_tid = falses(length(particle_tid))
#     for (tid, particle) in enumerate(particle_tid)
#         acceptedConstraints_tid[tid] = all(constraintsF(particle.simResults))
#     end
#     return accepted_tid .&& acceptedConstraints_tid
# end

# function testABCParticles(compareDataVSim::Function, particle_tid::Vector, dataMetrics, errorThresholds, constraintsF::Union{Nothing,Function}=nothing)
#     accepted_tid = falses(length(particle_tid))
#     for (tid,particle) in enumerate(particle_tid)
#         accepted_tid[tid] = testParticle(compareDataVSim, particle.simResults, dataMetrics, errorThresholds, constraintsF)
#     end
#     return accepted_tid
# end

# function testParticle(compareDataVSim::Function, simResults, dataMetrics, errorThresholds, constraintsF::Union{Nothing,Function}=nothing)
#     acceptedMetrics = compareDataVSim(dataMetrics, simResults, errorThresholds)
#     if isnothing(constraintsF)
#         return all(acceptedMetrics)
#     else
#         acceptedConstraints = constraintsF(simResults)
#         return all(acceptedMetrics) && all(acceptedConstraints)
#     end
# end

# function measureABCParticles(compareDataVSimError::Function, particle_tid::Vector, dataMetrics)
#     dists1 = compareDataVSimError(dataMetrics, particle_tid[1].simResults)
#     dists_tid = Vector{typeof(dists1)}(undef, length(particle_tid))
#     dists_tid[1] = dists1
#     for tid in 2:length(particle_tid)
#         dists_tid[tid] = compareDataVSimError(dataMetrics, particle_tid[tid].simResults)
#     end
#     return dists_tid
# end
