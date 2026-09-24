"""
    Particle(paramSet, simResults)

Package parameter values and existing simulation results into a particle.
No simulation is run; the supplied values are stored without copying.
`paramSet` may use any representation understood by the simulation, such as
a named tuple, tuple, vector, dictionary, or custom struct. Parameter indices
or keys should have consistent meanings across particles. Mutating a stored
parameter container also changes the values visible through the particle.
"""
struct Particle{P,R}
    paramSet::P
    simResults::R
end

"""
    packParticles(paramSets, simResults)

Package existing parameter sets and simulation results into particles.
Both vectors must have equal lengths. Entries are paired in iteration order,
and values are stored without copying, unless `referenceResults` and `referenceParams` are `false`.
Each element of `paramSets` holds one particle's parameters in any representation.
"""
function packParticles(
    paramSet_tid::AbstractVector,
    simResult_tid::AbstractVector;
    referenceResults::Bool=true,
    referenceParams::Bool=true,
)
    if length(paramSet_tid) != length(simResult_tid)
        throw(ArgumentError(
            "Expected equal numbers of parameter sets and simulation results; " *
            "got $(length(paramSet_tid)) and $(length(simResult_tid))."
        ))
    end

    if !referenceParams
        paramSet_tid = deepcopy(paramSet_tid)
    end
    if !referenceResults
        simResult_tid = deepcopy(simResult_tid)
    end

    return [Particle(params, result) for (params, result) in zip(paramSet_tid, simResult_tid)]
end

"""
    runABCParticles(
        runModelSim::Function,
        params_tid::AbstractVector,
        ctrlParams::Union{Tuple,Dict,NamedTuple}=NamedTuple(),
    )

Create multiple particles from a vector containing one parameter set per simulation.
Each parameter set is passed directly to `runModelSim`, which must understand
its representation. Parameter values are stored without copying.
"""
function runABCParticles(runModelSim::Function, params_tid::AbstractVector, ctrlParams::Union{Tuple,Dict,NamedTuple}=(;))
    return [
        runParticle(runModelSim, params, ctrlParams)
        for params in params_tid
    ]
end

"""
    runABCParticles(
        runModelSim::Function,
        params_tid_Pid::NamedTuple{Names, <:Tuple{Vararg{AbstractVector}}} where Names,
        ctrlParams::Union{Tuple,Dict,NamedTuple}=(;),
    )

Create multiple particles for the parameters in `params_tid_Pid`, which take the form of a `NamedTuple` of `Vector`s.
Each particle receives a named tuple with the same field names and the
corresponding entries from the parameter vectors, preserving their types.
"""
function runABCParticles(runModelSim::Function, params_tid_Pid::NamedTuple{Names, <:Tuple{Vararg{AbstractVector}}} where Names, ctrlParams::Union{Tuple,Dict,NamedTuple}=(;))
    nParticles = length(first(params_tid_Pid))
    for (pName, params_tid) in pairs(params_tid_Pid)
        if length(params_tid) != nParticles
            throw(ArgumentError("All parameter vectors in `params_tid_Pid` must have the same length. Parameter `$pName` has length $(length(params_tid)) while the first parameter has length $nParticles."))
        end
    end

    return [
        runParticle(runModelSim, map(values -> values[tid], params_tid_Pid), ctrlParams)
        for tid in 1:nParticles
    ]
end

function drawParams(priorDist_pid::Union{NamedTuple,Dict}, nParticles::Integer)
    pNames = keys(priorDist_pid) |> Tuple
    params_tid_Pid = NamedTuple{pNames}(Tuple(rand(dist, nParticles) for dist in values(priorDist_pid)))
    return params_tid_Pid
end

"""
    runABCParticles(
        runModelSim::Function,
        priorDist_pid::Union{NamedTuple,Dict},
        nParticles::Integer,
        ctrlParams::Union{Tuple,Dict,NamedTuple}=(;),
    )

Create multiple particles by first drawing `nParticles` parameters from the prior distributions in `priorDist_pid`.
"""
function runABCParticles(runModelSim::Function, priorDist_pid::Union{NamedTuple,Dict}, nParticles::Integer, ctrlParams::Union{Tuple,Dict,NamedTuple}=(;))
    params_tid_Pid = drawParams(priorDist_pid, nParticles)
    runABCParticles(runModelSim, params_tid_Pid, ctrlParams)
end

"""
    runParticle(runModelSim::Function, paramSet, ctrlParams)

Run `runModelSim(paramSet, ctrlParams)` and package the parameters and result
into a `Particle`. The simulation must understand the supplied parameter
representation. Parameters are passed and stored without copying.
"""
function runParticle(runModelSim::Function, pVal_pid, ctrlParams::Union{Tuple, Dict, NamedTuple})
    paramSet = pVal_pid # Use the passed parameter set directly
    simResults = runModelSim(paramSet, ctrlParams)
    @debug "Simulated particle" parameters=paramSet
    return Particle(
        paramSet,
        simResults
    )
end

# function runABC(runModelSim::Function, compareDataVSim::Function, priorDists_pid::NamedTuple, dataMetrics, ctrlParams::Dict, nParticles)

#     particle_tid = runABCParticles(runModelSim, priorDists_pid, nParticles, ctrlParams)

#     accepted_tid = testABCParticles(compareDataVSim, particle_tid, dataMetrics, ctrlParams[:thresholds])

#     return particle_tid, accepted_tid
# end

# function acceptedParams(particle_tid, accepted_tid)
#     [particle.paramSet for particle in particle_tid[accepted_tid]]
# end

# function rankedParams(particle_tid, tid_rank, number::Union{Int,Nothing}=nothing)
#     if isnothing(number)
#         return [particle.paramSet for particle in particle_tid[tid_rank]]
#     else
#         return [particle.paramSet for particle in particle_tid[tid_rank[1:number]]]
#     end
# end
