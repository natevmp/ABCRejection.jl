# ABCRejection

[![Build Status](https://github.com/natevmp/ABCRejection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/natevmp/ABCRejection.jl/actions/workflows/CI.yml?query=branch%3Amain)

Lightweight tools for performing [Approximate Bayesian Computation](https://en.wikipedia.org/wiki/Approximate_Bayesian_computation) in Julia.


## Usage

### 1. Generating Particles with `runABCParticles`

#### Define your Simulation Function
To generate simulations, first define a function to perform an instance of your model. The function must accept two arguments:

`runModelSim(paramSet, ctrlParams)`:
- `paramSet::NamedTuple`: The specific parameters sampled from your priors for a single simulation.
- `ctrlParams::Union{Dict, NamedTuple}`: Control parameters that remain constant across all simulations.

The function should return the result of the simulation that is later to be used for comparison with the reference data (e.g., a summary statistic or a time-series).

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
   Pass a `Vector` of `NamedTuple`s.
   ```julia
   params_pid_Tid = [(m = rand(), c = rand()) for _ in 1:100]
   particle_tid = runABCParticles(linearModel, params_pid_Tid, ctrlParams)
   ```

---

### 2. Ranking Particles with `rankParticles`

Instead of using a fixed error threshold, you can rank particles according to their distance from the observed data across multiple metrics.

#### Define your Distance Function
You must provide a function (e.g., `distDataVSim`) that calculates the distance. `distDataVSim(simResults, dataMetrics)`:
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

`rankParticles` sorts the particles in ascending order of their error (the first element is the "best" particle). If multiple metrics were used (as returned by the passed function `myDistance`), the particle acquires a rank for each metric. Its final ranking is the lowest of this set: $\mathrm{r} = \min \! \left( \left\lbrace \mathrm{r}_i \right\rbrace \right ) \forall i$.

```julia
dataMetrics = (5.0, 0.2) # Observed mean and std
tid_rank = rankParticles(myDistance, particle_tid, dataMetrics)
```
