module Bites

import Random, Distributions

function infect(human_prob::Float64, mosquito_prob::Float64, transmission_prob::Float64)::Bool
    Random.rand(Distributions.Bernoulli(human_prob*mosquito_prob*transmission_prob), 1)[1]
end

function one_way_bites(infecteds::Array{Float64, 1}, susceptibles::Array{Float64, 1}, transmission_prob::Float64)::Array{Bool, 1}
  result = zeros(Bool, length(susceptibles)) 
  if length(infecteds) > 0
    for i in 1:length(infecteds)
        susceptible_indexes = collect(1:length(susceptibles)) 
        [.!result]
        if length(susceptible_indexes) > 0
            for s in susceptible_indexes
                if infect(infecteds[i], susceptibles[s], transmission_prob)
                    result[s] = true
                end
            end
        end
    end
  end
  return result
end

export bite_steps
"""
    bite_steps(n_steps::Int64, n_humans::Int64, n_mosquitoes::Int64, human_infection_time::Int64, mosquito_life_span::Int64, human_probs::Array{Float64, 1}, mosquito_probs::Array{Float64, 1}, transmission_prob::Float64)::Tuple{Array{Int64, 1}, Array{Int64, 1}, Array{Int64, 1}}

Returns a Tuple containing (1) an array with number of infected mosquitoes at each time step, (2) an array with number of infected humans at each time step, and (3) an array with number of recovered humans at each time step.
"""
function bite_steps(n_steps::Int64, n_humans::Int64, n_mosquitoes::Int64, human_infection_time::Int64, mosquito_life_span::Int64, human_probs::Array{Float64, 1}, mosquito_probs::Array{Float64, 1}, transmission_prob::Float64)::Tuple{Array{Int64, 1}, Array{Int64, 1}, Array{Int64, 1}}

  # vectors of infections statusas follows: 0 = susceptible, >0 = infected, <0 = recovered. Everyone starts susceptible
  status_humans = zeros(Int8, n_humans)
  status_mosquitoes = zeros(Int8, n_mosquitoes)

  # vector of mosquito ages
  age_mosquitoes = Random.rand(Distributions.DiscreteUniform(0, mosquito_life_span), n_mosquitoes)

  # creating first infection
  status_humans[Random.rand(1:n_humans, 1)[1]] = 1

  n_human_infections = Vector{Int64}(undef, n_steps)
  n_mosquito_infections = Vector{Int64}(undef, n_steps)

  n_human_recovered = Vector{Int64}(undef, n_steps)

  n_human_infections[1] = 1
  n_mosquito_infections[1] = 0
  n_human_recovered[1] = 0


for s = 2:n_steps

    # update statuses based on time n_steps
    status_humans[findall(status_humans .> 0)] .+= 1
    status_humans[findall(status_humans .> human_infection_time)] .= -1

    # all mosquitoes age by 1 step
    age_mosquitoes .+= 1
    # mosquitoes over max age die and are reborn (under assumption of stable population) as susceptible
    status_mosquitoes[findall(age_mosquitoes .> mosquito_life_span)] .= 0
    age_mosquitoes[findall(age_mosquitoes .> mosquito_life_span)] .= 0

    # find indexes of infected humans (status>0)
    i_hs = findall(status_humans .> 0)
    # find indexes of susecptible mosquitoes (status == 0)
    s_ms = findall(status_mosquitoes .==0)
    
  # human-to-mosquito infections
    new_mosquito_infections = one_way_bites(human_probs[i_hs], mosquito_probs[s_ms], transmission_prob)
    if sum(new_mosquito_infections) > 0
        status_mosquitoes[s_ms[new_mosquito_infections]] .= 1
    end

    n_mosquito_infections[s] = sum(status_mosquitoes .> 0)

  s_hs = findall(status_humans .== 0)
  i_ms = findall(status_mosquitoes .==1)

    # mosquito-human infections
    new_human_infections = one_way_bites(mosquito_probs[i_ms], human_probs[s_hs], transmission_prob)
    if sum(new_human_infections) > 0
        status_humans[s_hs[new_human_infections]] .= 1
    end

    n_human_infections[s] = sum(status_humans.>0)
    n_human_recovered[s] = sum(status_humans.<0)
    
  end

  return n_mosquito_infections, n_human_infections, n_human_recovered

end

export distribute_bite_probabilities
"""
    distribute_bite_probabilities(D1, D2, N1::Int64, N2::Int64, expected_bites::Float64)::Tuple{Array{Float64, 1}, Array{Float64, 1}}

Returns a Tuple containing (1) per mosquito probabilities of being bitten for each person, and (2) per person probabilities of biting for each mosquito, with the per mosquito probabilities of being bitten for each person adjusted to the desired expected bites value.

## Parameters
* `D1` The distribution from which the human probabilities of being bitten should be drawn. A distribution from the Distributions.jl package. 
* `D2` The distribution from which the mosquito probabilities of biting should be drawn. A distribution from the Distributions.jl package. 
* `N1` The number of humans in the population. Int64.
* `N2` The number of mosquitoes in the population. Int64.
* `expected_bites` The expected number of total bites (total links in the bipartite network). Float64.
"""
function distribute_bite_probabilities(D1, D2, N1::Int64, N2::Int64, expected_bites::Float64)::Tuple{Array{Float64, 1}, Array{Float64, 1}}

  probs1 = Random.rand(D1, N1)
  probs2 = Random.rand(D2, N2)

  this_expected_bites = Float64(0)

  for i in probs1, j in probs2
    this_expected_bites += i*j
  end

  probs1 = probs1 .* (expected_bites/this_expected_bites)

  return probs1, probs2

end

export calculate_r0
"""
      calculate_r0(n_reps::Int64, infection_ts::Array{Int64,2}, seed_cases::Int64)

Returns a named Tuple containing (1) R0 (mean of reps), (2) R0 for each repetition, (3) mean R0 for successive numbers of repetitions.

## Parameters
* `n_reps` The number of repetitions in the simulation. In64. 
* `infection_ts` The matrix time series of number of infected cases at each time step. A two-dimensional array of In64 with dimensions equal to number of repetitions and number of time steps.
* `seed_cases` Number of infected cases seeded into the population at the beginning of each repetition. If `seed_cases=0`, then the first cases will be detected in the data. (For example, if we seed 1 infected human case into the population, then we would use `seed_cases=1` for the huiman R0 and `seed_cases=0` for the mosquito R0.)
"""
function calculate_r0(n_reps::Int64, infection_ts::Array{Int64,2}, seed_cases::Int64)

  R0s = Vector(undef, n_reps)

  if seed_cases > 0
    for i in 1:n_reps
      x = infection_ts[i,:]
      if maximum(x) > seed_cases
        R0s[i] = minimum(x[x.>seed_cases])
      else
        R0s[i] = 0
      end
    end  
  else
    for i in 1:n_reps
      x = n_mosquito_infections_reps[i,:]
      if maximum(x) > 1
        this_min = minimum(x[x.>0])
        if sum(x.>this_min)>0
          this_next_min = minimum(x[x.>this_min])
          R0s[i] = (this_next_min - this_min) / this_min
        else
          R0s[i] = 0
        end
      else
        R0s[i] = 0
      end
    end
  end

  R0 = mean(R0s)

  converge_check = Vector(undef, n_reps)
  for i in 1:n_reps
    converge_check[i] = mean(R0s[1:i])
  end

  return (R0=R0, R0_reps = R0s, converge_check = converge_check)

end



end
