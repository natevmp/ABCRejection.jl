using ABCRejection
using Test
using Distributions
using Random

mockModel(p, ctrl) = p.a + p.b

nParticles = 10
parName_pid = keys(priorDist_pid) |> Tuple
priorDist_pid = (a=Normal(1.,1.), b=Normal(2.,2.))
pVal_tid_pid = stack(rand(priorDist, nParticles) for priorDist in values(priorDist_pid))
params_tid = [NamedTuple{parName_pid}(Tuple(@view pVal_tid_pid[tid, :])) for tid in 1:nParticles]

particle_tid = runABCParticles(mockModel, priorDist_pid, 100)
particle_tid isa Vector{ABCRejection.Particle}
length(particle_tid)==100
typeof(particle_tid[1].paramSet)
particle_tid[1].simResults

for i in (a=1,b=2,c=3)
    println(i)
end

4!=3

length.(([1,2],[3,4],[5,6,7]))