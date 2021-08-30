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

end
