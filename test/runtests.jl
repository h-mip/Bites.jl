using Bites
using Test
import Random, Distributions

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

end
