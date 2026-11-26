using Bites
using Test
import Random, Distributions, Statistics

@testset "Bites.jl" begin
    Random.seed!(123)
    n_steps = 30
    n_reps = 10
    n_humans = 10
    n_mosquitoes = 40
    transmission_prob = .99
    human_infection_time = 3
    mosquito_life_span = 20
    mosquito_probs = fill(1/n_humans, n_mosquitoes)
    human_probs = fill(1/n_mosquitoes, n_humans)

    n_human_infections_reps = Array{Int64}(undef, n_reps, n_steps)
    n_mosquito_infections_reps = Array{Int64}(undef, n_reps, n_steps)
    n_human_recovered_reps = Array{Int64}(undef, n_reps, n_steps)
  
    # bite_steps
    for r in 1:n_reps
      n_mosquito_infections_reps[r, :], n_human_infections_reps[r, :], n_human_recovered_reps[r, :] = bite_steps(n_steps, n_humans, n_mosquitoes, human_infection_time, mosquito_life_span, human_probs, mosquito_probs, transmission_prob)
    end

    @test length(n_mosquito_infections_reps[1,:]) == n_steps
    @test length(n_human_infections_reps[1,:]) == n_steps
    @test length(n_human_recovered_reps[1,:]) == n_steps

    # distribute bites
    expected_bites = 2.0
    p1, p2 = distribute_bite_probabilities(Distributions.Uniform(0,1), Distributions.Uniform(0,1), 100, 400, expected_bites)
    this_expected_bites = Float64(0)
    for i in p1, j in p2
      this_expected_bites += i*j
    end
    @test abs(expected_bites - this_expected_bites) < .0001

    p1, p2 = distribute_bite_probabilities(Distributions.Normal(.9, 2), Distributions.Uniform(0,1), 100, 400, expected_bites)
    this_expected_bites = Float64(0)
    for i in p1, j in p2
      this_expected_bites += i*j
    end
    @test abs(expected_bites - this_expected_bites) < .0001

    p1, p2 = distribute_bite_probabilities(Distributions.Exponential(4), Distributions.Uniform(0,1), 100, 400, expected_bites)
    this_expected_bites = Float64(0)
    for i in p1, j in p2
      this_expected_bites += i*j
    end
    @test abs(expected_bites - this_expected_bites) < .0001

    p1, p2 = distribute_bite_probabilities(Distributions.Exponential(4), Distributions.Exponential(5), 100, 400, expected_bites)
    this_expected_bites = Float64(0)
    for i in p1, j in p2
      this_expected_bites += i*j
    end
    @test abs(expected_bites - this_expected_bites) < .0001

    # calculate_r0
    R0_result_h = calculate_r0(n_reps, n_human_infections_reps, 1)

    R0_result_m = calculate_r0(n_reps, n_mosquito_infections_reps, 0)
    
    @test length(R0_result_m.R0) == 1
    @test length(R0_result_m.R0_reps) == n_reps
    @test length(R0_result_m.converge_check) == n_reps
    @test length(R0_result_h.R0) == 1
    @test length(R0_result_h.R0_reps) == n_reps
    @test length(R0_result_h.converge_check) == n_reps

    # bite_steps_quad
    n_birds = 5
    n_horses = 2
    bird_infection_time = 3
    horse_infection_time = 4

    bird_probs = fill(0.8, n_birds)
    horse_probs = fill(0.5, n_horses)

    p_bird_to_mosquito = 0.9
    p_mosquito_to_bird = 0.9
    p_mosquito_to_human = 0.9
    p_mosquito_to_horse = 0.9

    Random.seed!(1234)
    n_mosq_q, n_bird_q, n_human_q, n_horse_q, n_bird_rec_q, n_human_rec_q, n_horse_rec_q = bite_steps_quad(
      n_steps, n_birds, n_mosquitoes, n_humans, n_horses, bird_infection_time, human_infection_time, horse_infection_time, mosquito_life_span,
      bird_probs, mosquito_probs, human_probs, horse_probs, p_bird_to_mosquito, p_mosquito_to_bird, p_mosquito_to_human, p_mosquito_to_horse;
      seed_birds=1, seed_mosquitoes=0, seed_humans=0, seed_horses=0)

    @test length(n_mosq_q) == n_steps
    @test length(n_bird_q) == n_steps
    @test length(n_human_q) == n_steps
    @test length(n_horse_q) == n_steps
    @test length(n_bird_rec_q) == n_steps
    @test length(n_human_rec_q) == n_steps
    @test length(n_horse_rec_q) == n_steps

    @test maximum(n_human_q) > 0
    @test maximum(n_mosq_q) > 0
    @test maximum(n_bird_q) > 0

    # dead-end check: seeding humans should not infect mosquitoes
    n_mosq_dead, n_bird_dead, _, _, _, _, _ = bite_steps_quad(
      n_steps, n_birds, n_mosquitoes, n_humans, n_horses, bird_infection_time, human_infection_time, horse_infection_time, mosquito_life_span,
      bird_probs, mosquito_probs, human_probs, horse_probs, p_bird_to_mosquito, p_mosquito_to_bird, p_mosquito_to_human, p_mosquito_to_horse;
      seed_birds=0, seed_mosquitoes=0, seed_humans=1, seed_horses=0)

    @test maximum(n_mosq_dead) == 0
    @test maximum(n_bird_dead) == 0

    @testset "Incubation delays" begin
      Random.seed!(321)
      n_steps = 6
      n_birds = 1
      n_mosquitoes = 1
      n_humans = 1
      n_horses = 0

      bird_probs = [1.0]
      mosquito_probs = [1.0]
      human_probs = [1.0]
      horse_probs = Float64[]

      p_bird_to_mosquito = 1.0
      p_mosquito_to_bird = 1.0
      p_mosquito_to_human = 1.0
      p_mosquito_to_horse = 1.0

      n_mosq_latent, n_bird_latent, n_human_latent, _, _, _, _ = bite_steps_quad(
        n_steps, n_birds, n_mosquitoes, n_humans, n_horses, 3, 3, 3, 10,
        bird_probs, mosquito_probs, human_probs, horse_probs,
        p_bird_to_mosquito, p_mosquito_to_bird, p_mosquito_to_human, p_mosquito_to_horse;
        seed_birds=1,
        seed_mosquitoes=0,
        bird_incubation_time=0,
        mosquito_incubation_time=2,
      )

      @test n_mosq_latent[2] == 1 # mosquito infected but latent
      @test n_human_latent[3] == 0 # cannot transmit while latent
      @test n_human_latent[4] >= 1 # becomes infectious after incubation
    end

    @testset "Human probability weighting" begin
      Random.seed!(11)
      n_steps = 4
      n_birds = 0
      n_mosquitoes = 1
      n_humans = 1
      n_horses = 0

      bird_probs = Float64[]
      mosquito_probs = [1.0]
      human_probs = [1.0]
      horse_probs = Float64[]

      p_bird_to_mosquito = 1.0
      p_mosquito_to_bird = 1.0
      p_mosquito_to_human = 1.0
      p_mosquito_to_horse = 1.0

      _, _, n_human_zero, _, _, _, _ = bite_steps_quad(
        n_steps, n_birds, n_mosquitoes, n_humans, n_horses, 1, 1, 1, 10,
        bird_probs, mosquito_probs, human_probs, horse_probs,
        p_bird_to_mosquito, p_mosquito_to_bird, p_mosquito_to_human, p_mosquito_to_horse;
        seed_birds=0,
        seed_mosquitoes=1,
        seed_humans=0,
        human_prob=0.0,
      )

      @test maximum(n_human_zero) == 0 # downweighted humans are not targeted

      _, _, n_human_default, _, _, _, _ = bite_steps_quad(
        n_steps, n_birds, n_mosquitoes, n_humans, n_horses, 1, 1, 1, 10,
        bird_probs, mosquito_probs, human_probs, horse_probs,
        p_bird_to_mosquito, p_mosquito_to_bird, p_mosquito_to_human, p_mosquito_to_horse;
        seed_birds=0,
        seed_mosquitoes=1,
        seed_humans=0,
      )

      @test maximum(n_human_default) > 0 # default behavior preserves human infections
    end

    @testset "bite_steps_quad_decay: independent host choice" begin
      Random.seed!(2026)
      nb = nh = no = 1000
      nm = 1000
      edges = BiteEdge[]
      counts = Dict{Int, Int}()
      out = bite_steps_quad_decay(40, nb, nm, nh, no, 7, 5, 5, 20,
        fill(0.5, nb), fill(1.0, nm), fill(0.5, nh), fill(0.1, no), 0.8, 0.8, 0.8, 0.8;
        seed_birds=1, bird_incubation_time=1, mosquito_incubation_time=3,
        gonotrophic_length=3, p_extra_bite=0.5, collect_bite_counts=true, bite_cycle_counts=counts, network_edges=edges)

      # positional destructuring of the first eight values still works
      ts_m, ts_b, ts_h, ts_o, rec_b, rec_h, rec_o, ever = out
      @test length(ts_b) == 40
      @test ts_b[1] == 1
      @test out.n_bird_infections === ts_b

      # bites land on host types in proportion to summed weights (0.5 : 0.5 : 0.1)
      share(t) = count(e -> e.target_type == t, edges) / length(edges)
      @test isapprox(share(:bird), 0.5 / 1.1; atol=0.02)
      @test isapprox(share(:human), 0.5 / 1.1; atol=0.02)
      @test isapprox(share(:horse), 0.1 / 1.1; atol=0.01)

      # feeding is spread over steps: a third of mosquitoes feed each step, two bites per meal on average
      per_step = [count(e -> e.step == s, edges) for s in 2:40]
      @test all(abs.(per_step .- nm / 3 * 2) .< 0.25 * nm / 3 * 2)

      # no host receives a pile of bites within one step
      worst = 0
      for s in 2:40
        load = Dict{Tuple{Symbol, Int}, Int}()
        for e in edges
          e.step == s || continue
          load[(e.target_type, e.target)] = get(load, (e.target_type, e.target), 0) + 1
        end
        worst = max(worst, maximum(values(load)))
      end
      @test worst <= 12

      # bites per meal are geometric with mean 1 / (1 - p_extra_bite), independent of the number of hosts
      meals = sum(values(counts))
      @test isapprox(sum(k * v for (k, v) in counts) / meals, 2.0; atol=0.1)
      @test isapprox(counts[1] / meals, 0.5; atol=0.05)
      @test !haskey(counts, 0)
      counts_small = Dict{Int, Int}()
      bite_steps_quad_decay(40, 100, nm, 50, 100, 7, 5, 5, 20,
        fill(0.5, 100), fill(1.0, nm), fill(0.5, 50), fill(0.1, 100), 0.8, 0.8, 0.8, 0.8;
        gonotrophic_length=3, p_extra_bite=0.5, collect_bite_counts=true, bite_cycle_counts=counts_small)
      @test isapprox(sum(k * v for (k, v) in counts_small) / sum(values(counts_small)), 2.0; atol=0.1)

      # generation bookkeeping adds up
      @test out.bird_generations[1] == 1
      @test sum(out.bird_generations) == ts_b[end] + rec_b[end]
      @test sum(out.mosquito_generations) == ever
      @test out.mosq_individuals_total >= nm
      @test isapprox(out.mosq_individuals_total, nm * (1 + 39 / 21); rtol=0.05)   # one cohort of 1000 turns over every 21 steps
      @test sum(out.human_generations) == ts_h[end] + rec_h[end]
      @test sum(out.horse_generations) == ts_o[end] + rec_o[end]

      # transmission log holds every non-seed infection
      tl = TransmissionEvent[]
      Random.seed!(7)
      out2 = bite_steps_quad_decay(60, 500, 500, 500, 500, 7, 5, 5, 20,
        fill(0.5, 500), fill(1.0, 500), fill(0.5, 500), fill(0.1, 500), 0.8, 0.8, 0.8, 0.8;
        bird_incubation_time=1, mosquito_incubation_time=3, gonotrophic_length=3, p_extra_bite=0.5, transmission_log=tl)
      @test length(tl) == (sum(out2.bird_generations) - 1) + sum(out2.mosquito_generations) + sum(out2.human_generations) + sum(out2.horse_generations)

      # per-individual weights: a zero-weight host is never bitten, a heavy one gets its share
      edges_w = BiteEdge[]
      wb = fill(1.0, 100)
      wb[1] = 0.0
      wb[2] = 50.0
      bite_steps_quad_decay(30, 100, 200, 0, 0, 7, 5, 5, 20, wb, fill(1.0, 200), Float64[], Float64[], 0.8, 0.8, 0.8, 0.8;
        gonotrophic_length=2, network_edges=edges_w)
      @test count(e -> e.target == 1, edges_w) == 0
      @test isapprox(count(e -> e.target == 2, edges_w) / length(edges_w), 50 / 148; atol=0.05)

      # dead-end hosts never start a chain
      out3 = bite_steps_quad_decay(30, 50, 100, 50, 50, 7, 5, 5, 20,
        fill(0.5, 50), fill(1.0, 100), fill(0.5, 50), fill(0.5, 50), 1.0, 1.0, 1.0, 1.0;
        seed_birds=0, seed_humans=1, seed_horses=1, gonotrophic_length=2)
      @test maximum(out3.n_mosquito_infections) == 0
      @test maximum(out3.n_bird_infections) == 0

      # removed / invalid keywords fail loudly
      args = (5, 10, 10, 10, 10, 7, 5, 5, 20, fill(0.5, 10), fill(1.0, 10), fill(0.5, 10), fill(0.1, 10), 0.8, 0.8, 0.8, 0.8)
      @test_throws ArgumentError bite_steps_quad_decay(args...; bite_decay=0.001)
      @test_throws ArgumentError bite_steps_quad_decay(args...; p_extra_bite=1.0)
      @test_throws ArgumentError bite_steps_quad_decay(args...; gonotrophic_length=0)
    end

    @testset "r0_analytic agrees with realised generations" begin
      N = 20_000
      M = 20_000
      ana = r0_analytic(; n_birds=N, n_mosquitoes=M, n_humans=N, n_horses=N,
        bird_weight=0.5, human_weight=0.5, horse_weight=0.1,
        p_bird_to_mosquito=0.8, p_mosquito_to_bird=0.8, p_mosquito_to_human=0.8, p_mosquito_to_horse=0.8,
        bird_infection_time=7, mosquito_life_span=20, mosquito_incubation_time=3, gonotrophic_length=3, p_extra_bite=0.5)
      @test isapprox(ana.share_bird, 0.5 / 1.1)
      @test isapprox(ana.bites_per_mosquito_per_step, 2 / 3)
      @test isapprox(ana.infectious_meals_per_infected_mosquito, 3.0)

      # large population so that depletion is negligible; 40 steps cover every first-generation infection
      nreps = 400
      g0_mosq = zeros(nreps)
      g1_birds = zeros(nreps)
      g1_humans = zeros(nreps)
      for r in 1:nreps
        Random.seed!(1000 + r)
        o = bite_steps_quad_decay(40, N, M, N, N, 7, 5, 5, 20,
          fill(0.5, N), fill(1.0, M), fill(0.5, N), fill(0.1, N), 0.8, 0.8, 0.8, 0.8;
          seed_birds=1, bird_incubation_time=1, mosquito_incubation_time=3, gonotrophic_length=3, p_extra_bite=0.5)
        g0_mosq[r] = length(o.mosquito_generations) >= 1 ? o.mosquito_generations[1] : 0
        g1_birds[r] = length(o.bird_generations) >= 2 ? o.bird_generations[2] : 0
        g1_humans[r] = length(o.human_generations) >= 2 ? o.human_generations[2] : 0
      end
      se(x) = Statistics.std(x) / sqrt(length(x))
      @test abs(Statistics.mean(g0_mosq) - ana.R_hv) < 4 * se(g0_mosq) + 0.05
      @test abs(Statistics.mean(g1_birds) - ana.R0) < 4 * se(g1_birds) + 0.1
      @test abs(Statistics.mean(g1_humans) - ana.spillover_human) < 4 * se(g1_humans) + 0.1
    end


    @testset "bite_steps_tracked: generations and seeding" begin
        Random.seed!(42)
        n_h, n_m = 200, 800
        hp, mp = distribute_bite_probabilities(Distributions.Exponential(1.0), 0.001:0.001, n_h, n_m, 300.0)
        res = bite_steps_tracked(200, n_h, n_m, 3, 20, hp, mp, 0.1)
        @test res.human_generations[1] == 1
        @test sum(res.human_generations) == res.n_human_infections[end] + res.n_human_recovered[end]
        @test length(res.n_human_infections) == 200
        @test_throws ArgumentError bite_steps_tracked(10, n_h, n_m, 3, 20, hp, mp, 0.1; seed=:nope)

        # degree-biased seeding picks better-connected index cases
        idx_u = [bite_steps_tracked(2, n_h, n_m, 3, 20, hp, mp, 0.1).index_human_prob for _ in 1:4000]
        idx_d = [bite_steps_tracked(2, n_h, n_m, 3, 20, hp, mp, 0.1; seed=:degree).index_human_prob for _ in 1:4000]
        @test isapprox(Statistics.mean(idx_u), Statistics.mean(hp); rtol=0.1)
        @test isapprox(Statistics.mean(idx_d), sum(abs2, hp) / sum(hp); rtol=0.1)

        # analytic expectation agrees with realised generation-1 counts (large, sparse population: little saturation)
        Random.seed!(7)
        n_h, n_m = 2000, 4000
        hp, mp = distribute_bite_probabilities(Distributions.Uniform(0.5, 1.5), 0.001:0.001, n_h, n_m, 4000.0)
        ana = r0_analytic_bipartite(hp, mp; transmission_prob=0.2, human_infection_time=3, mosquito_life_span=20)
        gen1(seed) = Statistics.mean(begin
            r = bite_steps_tracked(40, n_h, n_m, 3, 20, hp, mp, 0.2; seed=seed)
            length(r.human_generations) > 1 ? r.human_generations[2] : 0
        end for _ in 1:3000)
        @test isapprox(gen1(:uniform), ana.R0_uniform; rtol=0.12)
        @test isapprox(gen1(:degree), ana.R0_degree; rtol=0.12)
        @test ana.R0_degree > ana.R0_uniform
    end

end
