module ABCRejection

using Random, Distributions, Distances, Statistics

export runABCParticles, runParticle, testABCParticles, testParticle

struct Particle{T}
    paramSet::NamedTuple
    simResults::T
end

"""
    runABCParticles(
        runModelSim::Function,
        params_tid::Vector{<:NamedTuple},
        ctrlParams::Union{Dict,NamedTuple}=NamedTuple();
        verbose::Bool=false,
        randParams::Bool=true
    )

Create multiple particles for the parameters in `params_tid`
"""
function runABCParticles(runModelSim::Function, params_tid::Vector{<:NamedTuple}, ctrlParams::Union{Dict,NamedTuple}=NamedTuple(); verbose::Bool=false, randParams::Bool=false)
    particle_tid = Vector{Particle}(undef, length(params_tid))
    for tid in 1:length(params_tid)
        particle_tid[tid] = runParticle(runModelSim, params_tid[tid], ctrlParams; verbose, randParams)
    end
    return particle_tid
end

function runABCParticles(runModelSim::Function, priorDists_pid::NamedTuple, nParticles::Integer, ctrlParams::Dict; verbose::Bool=false, randParams::Bool=true)
    particle_tid = Vector{Particle}(undef, nParticles)
    for tid in 1:nParticles
        particle_tid[tid] = runParticle(runModelSim, priorDists_pid, ctrlParams; verbose, randParams)
    end
    return particle_tid
end

# function drawParams()
    
# end

function runParticle(runModelSim::Function, pDist_pid::NamedTuple, ctrlParams::Union{Dict, NamedTuple}; verbose::Bool=false, randParams::Bool=true)
    pVal_pid = 
        if randParams
            [rand(pDist) for pDist in pDist_pid]
        else
            pDist_pid
        end
    paramSet = (; zip(keys(pDist_pid), pVal_pid)...)
    simResults = runModelSim(paramSet, ctrlParams)
    if verbose
        println("parameter values of particle:")
        for (i,pid) in enumerate(keys(pDist_pid))
            println(string(pid)*": ", string(pVal_pid[i]))
        end
    end
    return Particle(
        paramSet,
        simResults
    )
end

function getParticleDistancesPerMetric(
        distDataVSim::Function,
        particle_tid::AbstractVector,
        dataMetrics;
    )
    nMetrics=length(dataMetrics)
    distance_tid_mid = Array{Float64,2}(undef, length(particle_tid), nMetrics)
    for (tid, particle) in enumerate(particle_tid)
        distance_tid_mid[tid, :] .= distDataVSim(particle.simResults, dataMetrics) |> collect
    end
    return distance_tid_mid
end

"""
    rankABCParticlesJointQuantiles(
        distDataVSim::Function,
        particle_tid::Vector,
        dataMetrics;
        nMetrics::Union{Nothing,Int}=nothing,
        verbose=false,
    )

Rank the particles according to their distance from the data in ascending order (first particle is closest).

`distDataVSim(particle.simResults, dataMetrics)`: Function to measure the distance between a single particle and the data. It must take two arguments, the first being the the `simResults` saved in each `particle` and the second the being all metrics of the data to which the particle is being compared.

"""
function rankABCParticlesJointQuantiles(
        distDataVSim::Function,
        particle_tid::Vector,
        dataMetrics::Tuple;
        # nMetrics::Union{Nothing,Int}=nothing,
        verbose=false,
    )
    # if isnothing(nMetrics) nMetrics=length(dataMetrics) end
    nMetrics=length(dataMetrics)
    distance_tid_mid = getParticleDistancesPerMetric(distDataVSim, particle_tid, dataMetrics)
    orderStat_tid_mid = Array{Int,2}(undef, (length(particle_tid),nMetrics))
    # for each metric, sort distances to get order statistic
    for mid in 1:nMetrics
        tid_tidSorted = sortperm(distance_tid_mid[:,mid])
        orderStat_tid_mid[:,mid] = invperm(tid_tidSorted)
    end
    # for each particle, take the maximum order statistic as new metric
    orderStatMax_tid = maximum(orderStat_tid_mid, dims=2)
    # sort particles by maximum order statistic
    tid_orderStatJoint = sortperm(orderStatMax_tid, dims=1) |> vec
    if verbose
        println("distances: distance_tid_mid = ") #! debug
        display(distance_tid_mid) #! debug
        println("order of distances: orderStat_tid_mid = ") #! debug
        display(orderStat_tid_mid) #! debug
        println("max distance per particle: orderStatMax_tid = ") #! debug
        display(orderStatMax_tid) #! debug
        println("order of Max rank: orderStatJoint_tid = ") #! debug
        display(tid_orderStatJoint) #! debug
    end
    return tid_orderStatJoint
end

function testABCParticlesQuantiles(
        distDataVSim::Function,
        particle_tid::Vector,
        dataMetrics,
        q::Real,
        constraintsF::Union{Nothing,Function}=nothing;
    )
    nMetrics=length(dataMetrics)
    distance_tid_mid = getParticleDistancesPerMetric(distDataVSim, particle_tid, dataMetrics)
    q_mid = [quantile((@view distance_tid_mid[:,mid]), q) for mid in 1:nMetrics]
    accepted_tid = [all(distance_tid_mid[tid,:].<q_mid) for tid in eachindex(particle_tid)]
    if isnothing(constraintsF)
        return accepted_tid
    end
    acceptedConstraints_tid = falses(length(particle_tid))
    for (tid, particle) in enumerate(particle_tid)
        acceptedConstraints_tid[tid] = all(constraintsF(particle.simResults))
    end
    return accepted_tid .&& acceptedConstraints_tid
end

function testABCParticles(compareDataVSim::Function, particle_tid::Vector, dataMetrics, errorThresholds, constraintsF::Union{Nothing,Function}=nothing)
    accepted_tid = falses(length(particle_tid))
    for (tid,particle) in enumerate(particle_tid)
        accepted_tid[tid] = testParticle(compareDataVSim, particle.simResults, dataMetrics, errorThresholds, constraintsF)
    end
    return accepted_tid
end

function testParticle(compareDataVSim::Function, simResults, dataMetrics, errorThresholds, constraintsF::Union{Nothing,Function}=nothing)
    acceptedMetrics = compareDataVSim(dataMetrics, simResults, errorThresholds)
    if isnothing(constraintsF)
        return all(acceptedMetrics)
    else
        acceptedConstraints = constraintsF(simResults)
        return all(acceptedMetrics) && all(acceptedConstraints)
    end
end

function measureABCParticles(compareDataVSimError::Function, particle_tid::Vector, dataMetrics)
    dists1 = compareDataVSimError(dataMetrics, particle_tid[1].simResults)
    dists_tid = Vector{typeof(dists1)}(undef, length(particle_tid))
    dists_tid[1] = dists1
    for tid in 2:length(particle_tid)
        dists_tid[tid] = compareDataVSimError(dataMetrics, particle_tid[tid].simResults)
    end
    return dists_tid
end

function runABC(runModelSim::Function, compareDataVSim::Function, priorDists_pid::NamedTuple, dataMetrics, ctrlParams::Dict, nParticles)

    particle_tid = runABCParticles(runModelSim, priorDists_pid, ctrlParams, nParticles)

    accepted_tid = testABCParticles(compareDataVSim, particle_tid, dataMetrics, ctrlParams[:thresholds])

    return particle_tid, accepted_tid
end

function acceptedParams(particle_tid, accepted_tid)
    [particle.paramSet for particle in particle_tid[accepted_tid]]
end

function rankedParams(particle_tid, tid_rank, number::Union{Int,Nothing}=nothing)
    if isnothing(number)
        return [particle.paramSet for particle in particle_tid[tid_rank]]
    else
        return [particle.paramSet for particle in particle_tid[tid_rank[1:number]]]
    end
end

end
