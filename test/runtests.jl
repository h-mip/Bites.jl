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

    # bite_steps
    n_mosquito_infections, n_human_infections, n_human_recovered = bite_steps(n_steps, n_humans, n_mosquitoes, human_infection_time, mosquito_life_span, human_probs, mosquito_probs, transmission_prob)
    @test length(n_mosquito_infections) == n_steps
    @test length(n_human_infections) == n_steps
    @test length(n_human_recovered) == n_steps

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

    # TO DO
    # calculate R0 using infection data generated above
    # R0_result_h = calculate_r0(n_steps, n_human_infections, seed_cases=1)

    # R0_result_m = calculate_r0(n_steps, n_mosquito_infections, seed_cases=0)

    @test 

end
