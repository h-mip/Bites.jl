using Bites
using Test
import Random, Distributions

@testset "Bites.jl" begin
    Random.seed!(123)
    n_steps = 3
    n_reps = 4
    n_humans = 10
    n_mosquitoes = 40
    transmission_prob = .99
    human_infection_time = 3
    mosquito_life_span = 20
    mosquito_probs = fill(1/n_humans, n_mosquitoes)
    human_probs = fill(1/n_mosquitoes, n_humans)
    n_mosquito_infections, n_human_infections, n_human_recovered = bite_steps(n_steps, n_humans, n_mosquitoes, human_infection_time, mosquito_life_span, human_probs, mosquito_probs, transmission_prob)
    @test length(n_mosquito_infections) == 3
    @test length(n_human_infections) == 3
    @test length(n_human_recovered) == 3

    expected_bites = 2.0
    p1, p2 = distribute_bite_probabilities(Distributions.Uniform(0,1), Distributions.Uniform(0,1), 100, 400, expected_bites)

    this_expected_bites = Float64(0)
    for i in p1, j in p2
      this_expected_bites += i*j
    end
  
    @test abs(expected_bites - this_expected_bites) < .0001

end
