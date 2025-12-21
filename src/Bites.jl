module Bites

import Random, Distributions, Statistics

function infect(human_prob::Float64, mosquito_prob::Float64, transmission_prob::Float64)::Bool
    this_p = human_prob * mosquito_prob * transmission_prob
    if this_p <= 0
      return false
    elseif this_p >= 1
      return true
    else
      return Random.rand(Distributions.Bernoulli(this_p))
    end
end

function one_way_bites(infecteds::Array{Float64, 1}, susceptibles::Array{Float64, 1}, transmission_prob::Float64)::Array{Bool, 1}
  result = zeros(Bool, length(susceptibles)) 
  if length(infecteds) > 0
    for i in 1:length(infecteds)
        susceptible_indexes = collect(1:length(susceptibles)) 
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

export bite_steps, bite_steps_quad, bite_steps_quad_decay
"""
    bite_steps(n_steps::Int64, n_humans::Int64, n_mosquitoes::Int64, human_infection_time::Int64, mosquito_life_span::Int64, human_probs::Array{Float64, 1}, mosquito_probs::Array{Float64, 1}, transmission_prob::Float64)::Tuple{Array{Int64, 1}, Array{Int64, 1}, Array{Int64, 1}}

Returns a Tuple containing (1) an array with number of infected mosquitoes at each time step, (2) an array with number of infected humans at each time step, and (3) an array with number of recovered humans at each time step.
"""
function bite_steps(n_steps::Int64, n_humans::Int64, n_mosquitoes::Int64, human_infection_time::Int64, mosquito_life_span::Int64, human_probs::Array{Float64, 1}, mosquito_probs::Array{Float64, 1}, transmission_prob::Float64)::Tuple{Array{Int64, 1}, Array{Int64, 1}, Array{Int64, 1}}

  # vectors of infections statusas follows: 0 = susceptible, >0 = infected, <0 = recovered. Everyone starts susceptible
  status_humans = zeros(Int, n_humans)
  status_mosquitoes = zeros(Int, n_mosquitoes)

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

"""
    bite_steps_quad(n_steps::Int64, n_birds::Int64, n_mosquitoes::Int64, n_humans::Int64, n_horses::Int64, bird_infection_time::Int64, human_infection_time::Int64, horse_infection_time::Int64, mosquito_life_span::Int64, bird_probs::Array{Float64, 1}, mosquito_probs::Array{Float64, 1}, human_probs::Array{Float64, 1}, horse_probs::Array{Float64, 1}, p_bird_to_mosquito::Float64, p_mosquito_to_bird::Float64, p_mosquito_to_human::Float64, p_mosquito_to_horse::Float64; seed_birds::Int=1, seed_mosquitoes::Int=0, seed_humans::Int=0, seed_horses::Int=0)

Simulates a quadripartite West Nile–style system with birds (reservoir), mosquitoes (vector), and dead-end hosts (humans and horses). Returns a tuple of infection and recovery time series ordered as (mosquito_infections, bird_infections, human_infections, horse_infections, bird_recovered, human_recovered, horse_recovered).
"""
function bite_steps_quad(n_steps::Int64, n_birds::Int64, n_mosquitoes::Int64, n_humans::Int64, n_horses::Int64, bird_infection_time::Int64, human_infection_time::Int64, horse_infection_time::Int64, mosquito_life_span::Int64, bird_probs::Array{Float64, 1}, mosquito_probs::Array{Float64, 1}, human_probs::Array{Float64, 1}, horse_probs::Array{Float64, 1}, p_bird_to_mosquito::Float64, p_mosquito_to_bird::Float64, p_mosquito_to_human::Float64, p_mosquito_to_horse::Float64; seed_birds::Int=1, seed_mosquitoes::Int=0, seed_humans::Int=0, seed_horses::Int=0)

  status_birds = zeros(Int, n_birds)
  status_mosquitoes = zeros(Int, n_mosquitoes)
  status_humans = zeros(Int, n_humans)
  status_horses = zeros(Int, n_horses)

  age_mosquitoes = Random.rand(Distributions.DiscreteUniform(0, mosquito_life_span), n_mosquitoes)

  if seed_birds > 0 && n_birds > 0
    status_birds[Random.rand(1:n_birds, min(seed_birds, n_birds))] .= 1
  end
  if seed_mosquitoes > 0 && n_mosquitoes > 0
    seeded_idx = Random.rand(1:n_mosquitoes, min(seed_mosquitoes, n_mosquitoes))
    status_mosquitoes[seeded_idx] .= 1
    for idx in seeded_idx
      if !infected_ever[idx]
        infected_ever[idx] = true
        mosq_ever_total += 1
      end
    end
  end
  if seed_humans > 0 && n_humans > 0
    status_humans[Random.rand(1:n_humans, min(seed_humans, n_humans))] .= 1
  end
  if seed_horses > 0 && n_horses > 0
    status_horses[Random.rand(1:n_horses, min(seed_horses, n_horses))] .= 1
  end


  n_bird_infections = Vector{Int64}(undef, n_steps)
  n_mosquito_infections = Vector{Int64}(undef, n_steps)
  n_human_infections = Vector{Int64}(undef, n_steps)
  n_horse_infections = Vector{Int64}(undef, n_steps)

  n_bird_recovered = Vector{Int64}(undef, n_steps)
  n_human_recovered = Vector{Int64}(undef, n_steps)
  n_horse_recovered = Vector{Int64}(undef, n_steps)

  n_bird_infections[1] = sum(status_birds .> 0)
  n_mosquito_infections[1] = sum(status_mosquitoes .> 0)
  n_human_infections[1] = sum(status_humans .> 0)
  n_horse_infections[1] = sum(status_horses .> 0)

  n_bird_recovered[1] = sum(status_birds .< 0)
  n_human_recovered[1] = sum(status_humans .< 0)
  n_horse_recovered[1] = sum(status_horses .< 0)

  for s = 2:n_steps

    status_birds[findall(status_birds .> 0)] .+= 1
    status_birds[findall(status_birds .> bird_infection_time)] .= -1

    status_humans[findall(status_humans .> 0)] .+= 1
    status_humans[findall(status_humans .> human_infection_time)] .= -1

    status_horses[findall(status_horses .> 0)] .+= 1
    status_horses[findall(status_horses .> horse_infection_time)] .= -1

    age_mosquitoes .+= 1
    status_mosquitoes[findall(age_mosquitoes .> mosquito_life_span)] .= 0
    age_mosquitoes[findall(age_mosquitoes .> mosquito_life_span)] .= 0

    i_bs = findall(status_birds .> 0)
    s_ms = findall(status_mosquitoes .== 0)
    if !isempty(i_bs) && !isempty(s_ms)
      new_mosquito_infections = one_way_bites(bird_probs[i_bs], mosquito_probs[s_ms], p_bird_to_mosquito)
      if sum(new_mosquito_infections) > 0
        status_mosquitoes[s_ms[new_mosquito_infections]] .= 1
      end
    end

    i_ms = findall(status_mosquitoes .> 0)

    s_bs = findall(status_birds .== 0)
    if !isempty(i_ms) && !isempty(s_bs)
      new_bird_infections = one_way_bites(mosquito_probs[i_ms], bird_probs[s_bs], p_mosquito_to_bird)
      if sum(new_bird_infections) > 0
        status_birds[s_bs[new_bird_infections]] .= 1
      end
    end

    s_hs = findall(status_humans .== 0)
    if !isempty(i_ms) && !isempty(s_hs)
      new_human_infections = one_way_bites(mosquito_probs[i_ms], human_probs[s_hs], p_mosquito_to_human)
      if sum(new_human_infections) > 0
        status_humans[s_hs[new_human_infections]] .= 1
      end
    end

    s_horses = findall(status_horses .== 0)
    if !isempty(i_ms) && !isempty(s_horses)
      new_horse_infections = one_way_bites(mosquito_probs[i_ms], horse_probs[s_horses], p_mosquito_to_horse)
      if sum(new_horse_infections) > 0
        status_horses[s_horses[new_horse_infections]] .= 1
      end
    end

    n_bird_infections[s] = sum(status_birds .> 0)
    n_mosquito_infections[s] = sum(status_mosquitoes .> 0)
    n_human_infections[s] = sum(status_humans .> 0)
    n_horse_infections[s] = sum(status_horses .> 0)

    n_bird_recovered[s] = sum(status_birds .< 0)
    n_human_recovered[s] = sum(status_humans .< 0)
    n_horse_recovered[s] = sum(status_horses .< 0)

  end

  return n_mosquito_infections, n_bird_infections, n_human_infections, n_horse_infections, n_bird_recovered, n_human_recovered, n_horse_recovered

end

"""
    bite_steps_quad_decay(
        n_steps::Int64,
        n_birds::Int64,
        n_mosquitoes::Int64,
        n_humans::Int64,
        n_horses::Int64,
        bird_infection_time::Int64,
        human_infection_time::Int64,
        horse_infection_time::Int64,
        mosquito_life_span::Int64,
        bird_probs::Array{Float64, 1},
        mosquito_probs::Array{Float64, 1},
        human_probs::Array{Float64, 1},
        horse_probs::Array{Float64, 1},
        p_bird_to_mosquito::Float64,
        p_mosquito_to_bird::Float64,
        p_mosquito_to_human::Float64,
        p_mosquito_to_horse::Float64;
        seed_birds::Int=1,
        seed_mosquitoes::Int=0,
        seed_humans::Int=0,
        seed_horses::Int=0,
        gonotrophic_length::Int=4,
        bite_decay::Float64=0.2,
        collect_bite_counts::Bool=false,
        bite_cycle_counts::Dict{Int, Int}=Dict{Int, Int}(),
        network_edges::Union{Nothing, Vector{NamedTuple{(:step, :mosquito, :target_type, :target, :success), Tuple{Int, Int, Symbol, Int, Bool}}}}=nothing,
    )

Optional logging: set `collect_bite_counts=true` and pass a `Dict{Int,Int}` via `bite_cycle_counts` to accumulate a histogram of bites per mosquito gonotrophic cycle. Provide `network_edges` as a vector of named tuples to capture bite edges for visualization; leave at `nothing` to avoid overhead. Each logged edge includes `success::Bool` to indicate whether transmission occurred. The function returns an eighth value: total mosquitoes ever infected (slot-level) over the run, useful for cumulative AR.
"""
function bite_steps_quad_decay(n_steps::Int64, n_birds::Int64, n_mosquitoes::Int64, n_humans::Int64, n_horses::Int64, bird_infection_time::Int64, human_infection_time::Int64, horse_infection_time::Int64, mosquito_life_span::Int64, bird_probs::Array{Float64, 1}, mosquito_probs::Array{Float64, 1}, human_probs::Array{Float64, 1}, horse_probs::Array{Float64, 1}, p_bird_to_mosquito::Float64, p_mosquito_to_bird::Float64, p_mosquito_to_human::Float64, p_mosquito_to_horse::Float64; seed_birds::Int=1, seed_mosquitoes::Int=0, seed_humans::Int=0, seed_horses::Int=0, gonotrophic_length::Int=4, bite_decay::Float64=0.2, collect_bite_counts::Bool=false, bite_cycle_counts::Dict{Int, Int}=Dict{Int, Int}(), network_edges::Union{Nothing, Vector{NamedTuple{(:step, :mosquito, :target_type, :target, :success), Tuple{Int, Int, Symbol, Int, Bool}}}}=nothing)

  status_birds = zeros(Int, n_birds)
  status_mosquitoes = zeros(Int, n_mosquitoes)
  status_humans = zeros(Int, n_humans)
  status_horses = zeros(Int, n_horses)

  age_mosquitoes = Random.rand(Distributions.DiscreteUniform(0, mosquito_life_span), n_mosquitoes)
  infected_ever = falses(n_mosquitoes)
  mosq_ever_total = 0

  if seed_birds > 0 && n_birds > 0
    status_birds[Random.rand(1:n_birds, min(seed_birds, n_birds))] .= 1
  end
  if seed_mosquitoes > 0 && n_mosquitoes > 0
    status_mosquitoes[Random.rand(1:n_mosquitoes, min(seed_mosquitoes, n_mosquitoes))] .= 1
  end
  if seed_humans > 0 && n_humans > 0
    status_humans[Random.rand(1:n_humans, min(seed_humans, n_humans))] .= 1
  end
  if seed_horses > 0 && n_horses > 0
    status_horses[Random.rand(1:n_horses, min(seed_horses, n_horses))] .= 1
  end

  record_bites = collect_bite_counts
  if record_bites
    empty!(bite_cycle_counts)
  end
  mosq_cycle_day = zeros(Int, n_mosquitoes)   # track days within gonotrophic cycle for all mosquitoes
  mosq_cycle_bites = zeros(Int, n_mosquitoes) # track bites within current cycle (used for decay even if not logging histogram)

  record_network = network_edges !== nothing
  if record_network
    empty!(network_edges)
  end

  n_bird_infections = Vector{Int64}(undef, n_steps)
  n_mosquito_infections = Vector{Int64}(undef, n_steps)
  n_human_infections = Vector{Int64}(undef, n_steps)
  n_horse_infections = Vector{Int64}(undef, n_steps)

  n_bird_recovered = Vector{Int64}(undef, n_steps)
  n_human_recovered = Vector{Int64}(undef, n_steps)
  n_horse_recovered = Vector{Int64}(undef, n_steps)

  n_bird_infections[1] = sum(status_birds .> 0)
  n_mosquito_infections[1] = sum(status_mosquitoes .> 0)
  n_human_infections[1] = sum(status_humans .> 0)
  n_horse_infections[1] = sum(status_horses .> 0)

  n_bird_recovered[1] = sum(status_birds .< 0)
  n_human_recovered[1] = sum(status_humans .< 0)
  n_horse_recovered[1] = sum(status_horses .< 0)

  cycle_len = max(gonotrophic_length, 1)
  bite_rate = 1 / cycle_len

  avg_m_prob = Statistics.mean(mosquito_probs)
  avg_b_prob = n_birds > 0 ? Statistics.mean(bird_probs) : 0.0
  avg_h_prob = n_humans > 0 ? Statistics.mean(human_probs) : 0.0
  avg_ho_prob = n_horses > 0 ? Statistics.mean(horse_probs) : 0.0
  scale_bird = (avg_m_prob * avg_b_prob) > 0 ? 1 / (avg_m_prob * avg_b_prob) : 0.0
  scale_human = (avg_m_prob * avg_h_prob) > 0 ? 1 / (avg_m_prob * avg_h_prob) : 0.0
  scale_horse = (avg_m_prob * avg_ho_prob) > 0 ? 1 / (avg_m_prob * avg_ho_prob) : 0.0

  for s = 2:n_steps

    mosq_cycle_day .+= 1
    completed = findall(>=(cycle_len), mosq_cycle_day)
    if !isempty(completed)
      if record_bites
        for idx in completed
          count = mosq_cycle_bites[idx]
          bite_cycle_counts[count] = get(bite_cycle_counts, count, 0) + 1
        end
      end
      mosq_cycle_bites[completed] .= 0
      mosq_cycle_day[completed] .= 0
    end

    status_birds[findall(status_birds .> 0)] .+= 1
    status_birds[findall(status_birds .> bird_infection_time)] .= -1

    status_humans[findall(status_humans .> 0)] .+= 1
    status_humans[findall(status_humans .> human_infection_time)] .= -1

    status_horses[findall(status_horses .> 0)] .+= 1
    status_horses[findall(status_horses .> horse_infection_time)] .= -1

    age_mosquitoes .+= 1
    dead_ms = findall(age_mosquitoes .> mosquito_life_span)
    if !isempty(dead_ms)
      if record_bites
        for idx in dead_ms
          count = mosq_cycle_bites[idx]
          bite_cycle_counts[count] = get(bite_cycle_counts, count, 0) + 1
        end
      end
      mosq_cycle_bites[dead_ms] .= 0
      mosq_cycle_day[dead_ms] .= 0
      status_mosquitoes[dead_ms] .= 0
      age_mosquitoes[dead_ms] .= 0
    end

    # Precompute shuffled orders to avoid ordering bias within the step
    bird_order = collect(1:n_birds); Random.shuffle!(bird_order)
    human_order = collect(1:n_humans); Random.shuffle!(human_order)
    horse_order = collect(1:n_horses); Random.shuffle!(horse_order)

    # Mosquito-centric biting with per-bite decay over the current gonotrophic cycle
    for m in 1:n_mosquitoes
      # Birds
      for b in bird_order
        if rand() < min(1.0, bite_rate * (bite_decay^mosq_cycle_bites[m]) * mosquito_probs[m] * bird_probs[b] * scale_bird)
          mosq_cycle_bites[m] += 1
          success = false
          if status_mosquitoes[m] > 0 && status_birds[b] == 0
            success = rand() < p_mosquito_to_bird
            if success
              status_birds[b] = 1
            end
          elseif status_birds[b] > 0 && status_mosquitoes[m] == 0
            success = rand() < p_bird_to_mosquito
            if success
              status_mosquitoes[m] = 1
              if !infected_ever[m]
                infected_ever[m] = true
                mosq_ever_total += 1
              end
            end
          end
          if record_network
            push!(network_edges, (step=s, mosquito=m, target_type=:bird, target=b, success=success))
          end
        end
      end

      # Humans
      for h in human_order
        if rand() < min(1.0, bite_rate * (bite_decay^mosq_cycle_bites[m]) * mosquito_probs[m] * human_probs[h] * scale_human)
          mosq_cycle_bites[m] += 1
          success = false
          if status_mosquitoes[m] > 0 && status_humans[h] == 0
            success = rand() < p_mosquito_to_human
            if success
              status_humans[h] = 1
            end
          end
          if record_network
            push!(network_edges, (step=s, mosquito=m, target_type=:human, target=h, success=success))
          end
        end
      end

      # Horses
      for ho in horse_order
        if rand() < min(1.0, bite_rate * (bite_decay^mosq_cycle_bites[m]) * mosquito_probs[m] * horse_probs[ho] * scale_horse)
          mosq_cycle_bites[m] += 1
          success = false
          if status_mosquitoes[m] > 0 && status_horses[ho] == 0
            success = rand() < p_mosquito_to_horse
            if success
              status_horses[ho] = 1
            end
          end
          if record_network
            push!(network_edges, (step=s, mosquito=m, target_type=:horse, target=ho, success=success))
          end
        end
      end
    end

    n_bird_infections[s] = sum(status_birds .> 0)
    n_mosquito_infections[s] = sum(status_mosquitoes .> 0)
    n_human_infections[s] = sum(status_humans .> 0)
    n_horse_infections[s] = sum(status_horses .> 0)

    n_bird_recovered[s] = sum(status_birds .< 0)
    n_human_recovered[s] = sum(status_humans .< 0)
    n_horse_recovered[s] = sum(status_horses .< 0)

  end

  return n_mosquito_infections, n_bird_infections, n_human_infections, n_horse_infections, n_bird_recovered, n_human_recovered, n_horse_recovered, mosq_ever_total

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
    this_expected_bites += i * j
  end

  if this_expected_bites == 0
    return probs1 .* 0, probs2 .* 0
  end

  scale = expected_bites / this_expected_bites
  probs1 = probs1 .* scale

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

  R0s = zeros(Float64, n_reps)

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
      x = infection_ts[i,:]
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

  R0 = Statistics.mean(R0s)

  converge_check = similar(R0s)
  for i in 1:n_reps
    converge_check[i] = Statistics.mean(R0s[1:i])
  end

  return (R0=R0, R0_reps = R0s, converge_check = converge_check)

end



end
