using ABCRejection
using Test
using Distributions
using Random

@testset "ABCRejection.jl" begin
    @testset "runParticle basic" begin

        mockModel(p, ctrl) = p.a + p.b
        inputParams = (a = 1.5, b = 2.5)
        
        particle = runParticle(mockModel, inputParams, (;))
        
        @test particle isa ABCRejection.Particle
        @test particle.simResults == 4.0
        @test particle.paramSet == inputParams
    end

    @testset "runABCParticles basics" begin
        mockModel(p, ctrl) = p.a + p.b

        priorDist_pid = (a=Normal(1.,1.), b=Normal(2.,2.))
        parName_pid = keys(priorDist_pid) |> Tuple
        nParticles = 100
        pVal_tid_pid = stack(rand(priorDist, nParticles) for priorDist in values(priorDist_pid))
        params_tid = [NamedTuple{parName_pid}(Tuple(@view pVal_tid_pid[tid, :])) for tid in 1:nParticles]

        particle_tid = runABCParticles(mockModel, priorDist_pid, nParticles)
        @test particle_tid isa Vector{ABCRejection.Particle}
        @test length(particle_tid)==nParticles
        @test particle_tid[1].paramSet isa @NamedTuple{a::Float64, b::Float64}
        @test particle_tid[1].simResults isa Float64

        particle_tid = runABCParticles(mockModel, params_tid)
        @test particle_tid isa Vector{ABCRejection.Particle}
        @test length(particle_tid)==nParticles
        @test particle_tid[1].paramSet isa @NamedTuple{a::Float64, b::Float64}
        @test particle_tid[1].simResults isa Float64

        pVal_tid_Pid = (a=rand(priorDist_pid.a, nParticles), b=rand(priorDist_pid.b, nParticles))
        particle_tid = runABCParticles(mockModel, pVal_tid_Pid)
        @test particle_tid isa Vector{ABCRejection.Particle}
        @test length(particle_tid)==nParticles
        @test particle_tid[1].paramSet isa @NamedTuple{a::Float64, b::Float64}
        @test particle_tid[1].simResults isa Float64
    end

end
