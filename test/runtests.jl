using ABCRejection
using Test
using Distributions
using Random

@testset "ABCRejection.jl" begin
    @testset "Ranking algorithm selection" begin
        # Metric ranks are (1,4,1), (2,2,2), (3,1,3), and (4,3,4).
        # Scaling the second metric distinguishes rank sums from distance sums.
        simResult_tid = [(1.0,400.0,1.0), (2.0,200.0,2.0), (3.0,100.0,3.0), (4.0,300.0,4.0)]
        particle_tid = [Particle((id=tid,), result) for (tid, result) in enumerate(simResult_tid)]
        distance(result, data) = result
        dataMetrics = (0.0, 0.0, 0.0)

        @test rankParticles(distance, particle_tid, dataMetrics) == [2, 3, 1, 4]
        @test rankParticles(distance, particle_tid, dataMetrics; algorithm=:max) == [2, 3, 1, 4]
        @test rankParticles(distance, particle_tid, dataMetrics; algorithm=:sum) == [1, 2, 3, 4]
        @test rankParticles(distance, particle_tid, (0.0,); nMetrics=3, algorithm=:sum) == [1, 2, 3, 4]

        tiedParticle_tid = [Particle((id=1,), (0.0,1.0)), Particle((id=2,), (0.0,0.0))]
        for algorithm in (:max, :sum)
            @test rankParticles(distance, tiedParticle_tid, (0.0,0.0); algorithm) == [2, 1]
            @test rankParticles((result, data) -> (1.0,), particle_tid, (0.0,); algorithm) == [1, 2, 3, 4]
            @test isempty(rankParticles(distance, empty(particle_tid), dataMetrics; algorithm))
        end

        @test_throws ArgumentError rankParticles(
            (result, data) -> error("Distance callback must not run"),
            particle_tid,
            dataMetrics;
            algorithm=:unknown,
        )
    end

    @testset "Maximum rank tie-breakers" begin
        # Resolve ties at the second and third worst ranks, respectively.
        @test ABCRejection.sortByMaxRank!([1 5 4; 2 3 5]) == [2, 1]
        @test ABCRejection.sortByMaxRank!([2 5 4; 1 4 5]) == [2, 1]
        # A better maximum takes priority over the sum or any later ranks.
        @test ABCRejection.sortByMaxRank!([4 4 4; 5 1 1]) == [1, 2]
        # Repeated worst ranks must be compared individually.
        @test ABCRejection.sortByMaxRank!([5 5 1; 5 4 4]) == [2, 1]
        # Metric labels do not distinguish identical sorted profiles.
        rank_tid_mid = [1 5 4; 4 1 5; 1 5 4]
        @test ABCRejection.sortByMaxRank!(rank_tid_mid) == [1, 2, 3]
        @test rank_tid_mid == [5 4 1; 5 4 1; 5 4 1]
        @test ABCRejection.sortByMaxRank!(reshape([3, 1, 1], :, 1)) == [2, 3, 1]
        @test ABCRejection.sortByMaxRank!([3 1 2]) == [1]
        @test isempty(ABCRejection.sortByMaxRank!(zeros(Int, 0, 3)))

        # Each column is a permutation of 1:5, so distances equal metric ranks.
        simResult_tid = [(2,5,4), (1,4,5), (3,1,1), (4,2,2), (5,3,3)]
        particle_tid = [Particle((id=tid,), result) for (tid, result) in enumerate(simResult_tid)]
        distance(result, data) = result
        dataMetrics = (0, 0, 0)
        @test rankParticles(distance, particle_tid, dataMetrics) == [3, 4, 5, 2, 1]
        @test rankParticles(distance, particle_tid, dataMetrics; algorithm=:max) == [3, 4, 5, 2, 1]
        @test rankParticles(distance, particle_tid, (0,); nMetrics=3) == [3, 4, 5, 2, 1]
        @test rankParticles(distance, particle_tid, dataMetrics; algorithm=:sum) == [3, 4, 2, 1, 5]
        @test rankParticles((result, data) -> (result[1], 100result[2], result[3]), particle_tid, dataMetrics) == [3, 4, 5, 2, 1]

        # Actual distance ties must receive shared ranks before aggregation.
        tiedParticle_tid = [Particle((id=1,), (1,2,3)), Particle((id=2,), (0,2,3))]
        @test rankParticles(distance, tiedParticle_tid, dataMetrics) == [2, 1]
        symmetricParticle_tid = [Particle((id=1,), (1,2)), Particle((id=2,), (2,1))]
        @test rankParticles(distance, symmetricParticle_tid, (0,0)) == [1, 2]
    end

    @testset "Shared competition rank buffers" begin
        distance_tid_mid = [4.0 2.0; 1.0 5.0; 1.0 2.0; 9.0 1.0]
        originalDistance_tid_mid = copy(distance_tid_mid)
        rank_tid = Vector{Int}(undef, 4)
        tid_tidSorted = Vector{Int}(undef, 4)

        @test ABCRejection.competitionRanks!(rank_tid, @view(distance_tid_mid[:, 1]), tid_tidSorted) === rank_tid
        @test rank_tid == [3, 1, 1, 4]
        ABCRejection.competitionRanks!(rank_tid, @view(distance_tid_mid[:, 2]), tid_tidSorted)
        @test rank_tid == [2, 4, 2, 1]
        @test distance_tid_mid == originalDistance_tid_mid

        @test ABCRejection.rankBySum(distance_tid_mid) == [3, 1, 2, 4]
        @test ABCRejection.rankByMax(distance_tid_mid) == [3, 1, 2, 4]
        @test distance_tid_mid == originalDistance_tid_mid
    end

    @testset "Ranking helpers against direct-count reference" begin
        rng = MersenneTwister(42)
        for nParticles in (0, 1, 12), nMetrics in (1, 3, 8)
            distance_tid_mid = rand(rng, 0:4, nParticles, nMetrics)
            # Independent reference: count strictly smaller values per metric.
            referenceRank_tid_mid = [
                1 + count(value -> value < distance_tid_mid[tid, mid], @view(distance_tid_mid[:, mid]))
                for tid in 1:nParticles, mid in 1:nMetrics
            ]
            expectedSum_tid = vec(sum(referenceRank_tid_mid; dims=2))
            expectedProfile_tid = [Tuple(sort(collect(row); rev=true)) for row in eachrow(referenceRank_tid_mid)]
            @test ABCRejection.rankBySum(distance_tid_mid) == sortperm(expectedSum_tid)
            @test ABCRejection.rankByMax(distance_tid_mid) == sortperm(expectedProfile_tid)
        end
    end

    @testset "packParticles" begin
        paramSets = [(a=1.0, b=2.0), (a=3.0, b=4.0)]
        simResults = [[0.1, 0.2], [0.3, 0.4]]

        particles = packParticles(paramSets, simResults)
        @test length(particles) == 2
        @test all(particle isa Particle for particle in particles)
        @test [particle.paramSet for particle in particles] == paramSets
        @test all(particles[i].simResults === simResults[i] for i in eachindex(particles))
        @test isempty(packParticles(empty(paramSets), empty(simResults)))
        @test_throws ArgumentError packParticles(paramSets, simResults[1:1])
        @test_throws ArgumentError packParticles(paramSets[1:1], simResults)
    end

    @testset "Positional parameter containers" begin
        positionalModel(p, ctrl) = p[1] + p[2] + ctrl.offset
        distance(result, data) = (abs(result - data[1]),)
        ctrlParams = (offset=0.5,)

        for paramSets in ([(1.0, 2), (3.0, 4)], [[1.0, 2.0], [3.0, 4.0]])
            simResults = [3.5, 7.5]
            particle = Particle(paramSets[1], simResults[1])
            @test particle isa Particle{eltype(paramSets),Float64}
            @test particle.paramSet === paramSets[1]

            # Views exercise AbstractVector support for the outer collections.
            particles = packParticles(view(paramSets, :), view(simResults, :))
            @test eltype(particles) == typeof(particle)
            @test [p.simResults for p in particles] == simResults
            @test all(particles[i].paramSet === paramSets[i] for i in eachindex(paramSets))
            @test_throws ArgumentError packParticles(paramSets, simResults[1:1])

            simulated = runParticle(positionalModel, paramSets[1], ctrlParams)
            @test simulated.paramSet === paramSets[1]
            @test simulated.simResults == simResults[1]

            simulatedParticles = runABCParticles(positionalModel, view(paramSets, :), ctrlParams)
            @test eltype(simulatedParticles) == typeof(particle)
            @test [p.simResults for p in simulatedParticles] == simResults
            @test all(simulatedParticles[i].paramSet === paramSets[i] for i in eachindex(paramSets))
            @test rankParticles(distance, simulatedParticles, (6.0,)) == [2, 1]
        end
    end

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
        @test particle_tid isa Vector{<:ABCRejection.Particle}
        @test isconcretetype(eltype(particle_tid))
        @test length(particle_tid)==nParticles
        @test particle_tid[1].paramSet isa @NamedTuple{a::Float64, b::Float64}
        @test particle_tid[1].simResults isa Float64

        particle_tid = runABCParticles(mockModel, params_tid)
        @test particle_tid isa Vector{<:ABCRejection.Particle}
        @test length(particle_tid)==nParticles
        @test particle_tid[1].paramSet isa @NamedTuple{a::Float64, b::Float64}
        @test particle_tid[1].simResults isa Float64

        pVal_tid_Pid = (a=rand(priorDist_pid.a, nParticles), b=rand(priorDist_pid.b, nParticles))
        particle_tid = runABCParticles(mockModel, pVal_tid_Pid)
        @test particle_tid isa Vector{<:ABCRejection.Particle}
        @test length(particle_tid)==nParticles
        @test particle_tid[1].paramSet isa @NamedTuple{a::Float64, b::Float64}
        @test particle_tid[1].simResults isa Float64
    end

    @testset "runABCParticles named parameter columns" begin
        params_tid_Pid = (b=[2, 4], a=[1.5, 3.5])
        expectedParams_tid = [(b=2, a=1.5), (b=4, a=3.5)]
        nCalls = Ref(0)
        function columnModel(params, ctrl)
            nCalls[] += 1
            return params.a + params.b + ctrl.offset
        end
        ctrlParams = (offset=0.5,)

        particle_tid = runABCParticles(columnModel, params_tid_Pid, ctrlParams)
        @test nCalls[] == 2
        @test eltype(particle_tid) == Particle{eltype(expectedParams_tid),Float64}
        @test [p.paramSet for p in particle_tid] == expectedParams_tid
        @test [p.simResults for p in particle_tid] == [4.0, 8.0]

        emptyColumns = (b=Int[], a=Float64[])
        @test isempty(runABCParticles(columnModel, emptyColumns, ctrlParams))
        @test nCalls[] == 2
        @test_throws ArgumentError runABCParticles(columnModel, (b=[2], a=[1.5, 3.5]), ctrlParams)
        @test nCalls[] == 2

        # Column entries may themselves be containers; preserve their references.
        nestedColumns = (a=[[1.0, 2.0], [3.0, 4.0]], b=[2, 4])
        nestedParticle_tid = runABCParticles((p, ctrl) -> sum(p.a) + p.b, nestedColumns)
        @test isconcretetype(eltype(nestedParticle_tid))
        @test all(nestedParticle_tid[i].paramSet.a === nestedColumns.a[i] for i in 1:2)
        @test [p.simResults for p in nestedParticle_tid] == [5.0, 11.0]
    end

    @testset "runABCParticles preserves complete parameter sets" begin
        for paramSet_tid in (
            [(a=1.0, b=2), (a=3.0, b=4)],
            [(1.0, 2), (3.0, 4)],
            [[1.0, 2.0], [3.0, 4.0]],
            [Dict(:a=>1.0, :b=>2.0), Dict(:a=>3.0, :b=>4.0)],
        )
            particle_tid = runABCParticles((params, ctrl) -> params, paramSet_tid)
            @test length(particle_tid) == length(paramSet_tid)
            @test isconcretetype(eltype(particle_tid))
            @test all(particle_tid[i].paramSet === paramSet_tid[i] for i in eachindex(paramSet_tid))
            @test all(particle_tid[i].simResults === paramSet_tid[i] for i in eachindex(paramSet_tid))
        end
    end

end
