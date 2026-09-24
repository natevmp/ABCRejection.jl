module ABCRejection

using Random: rand

export Particle, packParticles, runABCParticles, runParticle, rankParticles

include("particles.jl")
include("ranking.jl")

end
