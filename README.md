# ABCRejection

[![Build Status](https://github.com/natevmp/ABCRejection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/natevmp/ABCRejection.jl/actions/workflows/CI.yml?query=branch%3Amain)

Lightweight tools for performing [Approximate Bayesian Computation](https://en.wikipedia.org/wiki/Approximate_Bayesian_computation) in Julia.


## Usage

### 1. Generating Particles with `runABCParticles`

#### Define your Simulation Function
To generate simulations, first define a function to perform an instance of your model. The function must accept two arguments:

`runModelSim(paramSet, ctrlParams)`:
- `paramSet`: Parameters for a single simulation, in a representation your model understands, such as a named tuple, tuple, vector, dictionary, or custom struct.
- `ctrlParams::Union{Dict, NamedTuple}`: Control parameters that remain constant across all simulations.

The function should return the result of the simulation that is later to be used for comparison with the reference data (e.g., a summary statistic or a time-series).

Parameter indices or keys should identify the same parameters across particles. Parameters are passed and stored without copying, so later mutations to a parameter container are visible through its particle. Prior sampling and the named-tuple-of-vectors input below produce named-tuple parameter sets.

```julia
function linearModel(paramSet, ctrlParams)
    # Example: a simple linear model y = mx + c
    # paramSet contains m and c, ctrlParams contains x
    result = paramSet.m .* ctrlParams.x .+ paramSet.c
    return result
end
```

#### Running the Simulations
There are three ways to call `runABCParticles`, depending on how you to handle your parameters:

1. **From Priors (Automatic Drawing):**
   Pass your prior distributions and the desired number of particles. The package will handle the sampling for you.
   ```julia
   priorDist_pid = (m = Normal(0, 1), c = Normal(10, 5))
   particle_tid = runABCParticles(linearModel, priorDist_pid, 100, ctrlParams)
   ```

2. **Pre-drawn Parameters (Column-major):**
   Pass a `NamedTuple` where each field is a vector of pre-drawn values.
   ```julia
   params_tid_Pid = (m = rand(100), c = rand(100))
   particle_tid = runABCParticles(linearModel, params_tid_Pid, ctrlParams)
   ```

3. **Pre-drawn Parameters (Row-major):**
   Pass an `AbstractVector` containing one parameter set per simulation. Each set is passed directly to your model. For the named-tuple model above:
   ```julia
   params_pid_Tid = [(m = rand(), c = rand()) for _ in 1:100]
   particle_tid = runABCParticles(linearModel, params_pid_Tid, ctrlParams)
   ```

   A model using positional indexing can instead accept tuples or vectors:
   ```julia
   positionalModel(params, ctrl) = params[1] .* ctrl.x .+ params[2]
   params_pid_Tid = [(1.0, 2.0), (3.0, 4.0)] # or [[1.0, 2.0], [3.0, 4.0]]
   particle_tid = runABCParticles(positionalModel, params_pid_Tid, ctrlParams)
   ```

#### Packaging Existing Results

Use `Particle(paramSet, simResults)` for one existing simulation, or `packParticles(paramSets, simResults)` for vectors of parameter sets and results. The vectors must have equal lengths; entries are paired in iteration order and stored without copying.

---

### 2. Ranking Particles with `rankParticles`

Instead of using a fixed error threshold, you can rank particles according to their distance from the observed data across multiple metrics.

#### Define your Distance Function
You must provide a function (e.g., `distanceFunc`) that calculates the distance. `distanceFunc(simResults, dataMetrics)`:
- `simResults`: The output returned by your `runModelSim` function for a single particle.
- `dataMetrics`: The observed data or target metrics you are comparing against. If `length(dataMetrics)>1`, (e.g. as a `Tuple` or `Vector`), the data will be compared to multiple metrics, where the number of metrics is equal to `length(dataMetrics)`.



The function must return a distance (or a vector of distances if using multiple metrics) between the simulation and the data.
```julia
function myDistance(simResults, dataMetrics)
    # Returns a tuple of distances for two different metrics
    dist1 = abs(mean(simResults) - dataMetrics[1])
    dist2 = abs(std(simResults) - dataMetrics[2])
    return (dist1, dist2)
end
```

`rankParticles` assigns competition ranks separately for each metric: smaller distances receive smaller ranks, and equal distances share a rank. The `algorithm` keyword determines how those ranks are combined:

- `algorithm=:max` (default): compare each particle's ranks from worst to best. Compare the maximum rank first; if tied, compare the second-largest rank, then the third-largest, and so on.
- `algorithm=:sum`: use the sum of the ranks across metrics.

Lower ranks or sums are better. For `:max`, particles with identical sorted rank profiles retain input order. For `:sum`, equal sums retain input order without a secondary tie-breaker. Both algorithms return particle indices rather than changing the particle vector.

Signed zeros (`-0.0` and `0.0`) share a rank, and infinite distances are allowed.
Particles with `NaN` in any distance are excluded before ranks are computed,
with one warning listing their positions. The returned indices refer to the
supplied particle collection and omit excluded particles; if all particles are
excluded, the result is empty.

```julia
dataMetrics = (5.0, 0.2) # Observed mean and std
tid_rank = rankParticles(myDistance, particle_tid, dataMetrics)
tid_rankSum = rankParticles(myDistance, particle_tid, dataMetrics; algorithm=:sum)
```
