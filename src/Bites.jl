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
    status_mosquitoes[Random.rand(1:n_mosquitoes, min(seed_mosquitoes, n_mosquitoes))] .= 1
  end
  if seed_humans > 0 && n_humans > 0
    status_humans[Random.rand(1:n_humans, min(seed_humans, n_humans))] .= 1
  end
  if seed_horses > 0 && n_horses > 0
    status_horses[Random.rand(1:n_horses, min(seed_horses, n_horses))] .= 1
  end

  # Optional per-cycle bite counting (aggregated histogram only; avoids large allocations)
  record_bites = collect_bite_counts
  if record_bites
    empty!(bite_cycle_counts)
    mosq_cycle_day = zeros(Int, n_mosquitoes)
    mosq_cycle_bites = zeros(Int, n_mosquitoes)
  else
    mosq_cycle_day = nothing
    mosq_cycle_bites = nothing
  end

  # Optional bite-level edge logging for network visualization; caller can toggle for a single run
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

  if record_bites
    for count in mosq_cycle_bites
      bite_cycle_counts[count] = get(bite_cycle_counts, count, 0) + 1
    end
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
        network_edges::Union{Nothing, Vector{NamedTuple{(:step, :mosquito, :target_type, :target), Tuple{Int, Int, Symbol, Int}}}}=nothing,
    )

Quadripartite simulation where each infected mosquito's chance to take additional bites within a day decays after each bite. A shorter `gonotrophic_length` increases bite attempts per day. The baseline functions remain unchanged; use this variant when you want per-day bite tapering.

Optional logging: set `collect_bite_counts=true` and pass a `Dict{Int,Int}` via `bite_cycle_counts` to accumulate a histogram of bites per mosquito gonotrophic cycle. Provide `network_edges` as a vector of named tuples to capture bite edges for visualization; leave at `nothing` to avoid overhead.
"""
function bite_steps_quad_decay(n_steps::Int64, n_birds::Int64, n_mosquitoes::Int64, n_humans::Int64, n_horses::Int64, bird_infection_time::Int64, human_infection_time::Int64, horse_infection_time::Int64, mosquito_life_span::Int64, bird_probs::Array{Float64, 1}, mosquito_probs::Array{Float64, 1}, human_probs::Array{Float64, 1}, horse_probs::Array{Float64, 1}, p_bird_to_mosquito::Float64, p_mosquito_to_bird::Float64, p_mosquito_to_human::Float64, p_mosquito_to_horse::Float64; seed_birds::Int=1, seed_mosquitoes::Int=0, seed_humans::Int=0, seed_horses::Int=0, gonotrophic_length::Int=4, bite_decay::Float64=0.2, collect_bite_counts::Bool=false, bite_cycle_counts::Dict{Int, Int}=Dict{Int, Int}(), network_edges::Union{Nothing, Vector{NamedTuple{(:step, :mosquito, :target_type, :target), Tuple{Int, Int, Symbol, Int}}}}=nothing)

  status_birds = zeros(Int, n_birds)
  status_mosquitoes = zeros(Int, n_mosquitoes)
  status_humans = zeros(Int, n_humans)
  status_horses = zeros(Int, n_horses)

  age_mosquitoes = Random.rand(Distributions.DiscreteUniform(0, mosquito_life_span), n_mosquitoes)

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

  for s = 2:n_steps

    if record_bites
      mosq_cycle_day .+= 1
      completed = findall(>=(cycle_len), mosq_cycle_day)
      if !isempty(completed)
        for idx in completed
          count = mosq_cycle_bites[idx]
          bite_cycle_counts[count] = get(bite_cycle_counts, count, 0) + 1
        end
        mosq_cycle_bites[completed] .= 0
        mosq_cycle_day[completed] .= 0
      end
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
        mosq_cycle_bites[dead_ms] .= 0
        mosq_cycle_day[dead_ms] .= 0
      end
      status_mosquitoes[dead_ms] .= 0
      age_mosquitoes[dead_ms] .= 0
    end

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
    s_hs = findall(status_humans .== 0)
    s_horses = findall(status_horses .== 0)

    # mosquito -> bird/human/horse with per-bite decay
    for m in i_ms
      attempt = 0
      while rand() < (bite_rate * (bite_decay^attempt))
        if record_bites
          mosq_cycle_bites[m] += 1
        end
        w_b = isempty(s_bs) ? 0.0 : sum(bird_probs[s_bs])
        w_h = isempty(s_hs) ? 0.0 : sum(human_probs[s_hs])
        w_ho = isempty(s_horses) ? 0.0 : sum(horse_probs[s_horses])
        total_w = w_b + w_h + w_ho
        total_w == 0 && break

        r = rand() * total_w
        if r < w_b
          # pick bird
          if length(s_bs) == 1
            idx = s_bs[1]
          else
            target = rand() * w_b
            accum = 0.0
            idx = s_bs[1]
            for j in s_bs
              accum += bird_probs[j]
              if accum >= target
                idx = j
                break
              end
            end
          end
          if record_network
            push!(network_edges, (step=s, mosquito=m, target_type=:bird, target=idx))
          end
          if rand() < p_mosquito_to_bird
            status_birds[idx] = 1
            deleteat!(s_bs, findfirst(==(idx), s_bs))
          end
        elseif r < w_b + w_h
          # pick human
          if length(s_hs) == 1
            idx = s_hs[1]
          else
            target = rand() * w_h
            accum = 0.0
            idx = s_hs[1]
            for j in s_hs
              accum += human_probs[j]
              if accum >= target
                idx = j
                break
              end
            end
          end
          if record_network
            push!(network_edges, (step=s, mosquito=m, target_type=:human, target=idx))
          end
          if rand() < p_mosquito_to_human
            status_humans[idx] = 1
            deleteat!(s_hs, findfirst(==(idx), s_hs))
          end
        else
          # pick horse
          if length(s_horses) == 1
            idx = s_horses[1]
          else
            target = rand() * w_ho
            accum = 0.0
            idx = s_horses[1]
            for j in s_horses
              accum += horse_probs[j]
              if accum >= target
                idx = j
                break
              end
            end
          end
          if record_network
            push!(network_edges, (step=s, mosquito=m, target_type=:horse, target=idx))
          end
          if rand() < p_mosquito_to_horse
            status_horses[idx] = 1
            deleteat!(s_horses, findfirst(==(idx), s_horses))
          end
        end

        attempt += 1
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
