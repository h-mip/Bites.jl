module Bites

import Random, Distributions, Statistics

"""
    infect(human_prob::Float64, mosquito_prob::Float64, transmission_prob::Float64)::Bool

Determines whether a transmission event occurs given the probability of an infected host being bitten, the probability of a susceptible host being bitten, and the transmission probability.

Returns `true` if transmission occurs, `false` otherwise.

## Parameters
* `human_prob` - Probability that the human (or other host) is involved in a bite event. Float64.
* `mosquito_prob` - Probability that the mosquito is involved in a bite event. Float64.
* `transmission_prob` - Probability of transmission given a bite event between infected and susceptible individuals. Float64.
"""
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

"""
    one_way_bites(infecteds::Array{Float64, 1}, susceptibles::Array{Float64, 1}, transmission_prob::Float64)::Array{Bool, 1}

Simulates one-way transmission from infected individuals to susceptible individuals.

Returns a boolean array indicating which susceptible individuals become infected.

## Parameters
* `infecteds` - Array of probabilities for infected individuals. Float64 array.
* `susceptibles` - Array of probabilities for susceptible individuals. Float64 array.
* `transmission_prob` - Probability of transmission given contact. Float64.
"""
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

  # vectors of infections statuses as follows: 0 = susceptible, >0 = infected, <0 = recovered. Everyone starts susceptible
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
    # find indexes of susceptible mosquitoes (status == 0)
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
  bite_steps_quad(
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
    bird_incubation_time::Int=10,
    mosquito_incubation_time::Int=4,
    human_prob::Union{Float64, AbstractVector{<:Float64}, Nothing}=nothing,
  )

Simulates a quadripartite West Nile–style system with birds (reservoir), mosquitoes (vector), and dead-end hosts (humans and horses). Newly infected birds and mosquitoes enter a latent state for `bird_incubation_time` or `mosquito_incubation_time` steps before becoming infectious; infectious duration then follows the respective `*_infection_time`. Optionally weight humans via `human_prob` (scalar or per-human vector) when normalizing mosquito target selection across bird/human/horse hosts. Returns a tuple of infection and recovery time series ordered as (mosquito_infections, bird_infections, human_infections, horse_infections, bird_recovered, human_recovered, horse_recovered).
"""
function bite_steps_quad(n_steps::Int64, n_birds::Int64, n_mosquitoes::Int64, n_humans::Int64, n_horses::Int64, bird_infection_time::Int64, human_infection_time::Int64, horse_infection_time::Int64, mosquito_life_span::Int64, bird_probs::Array{Float64, 1}, mosquito_probs::Array{Float64, 1}, human_probs::Array{Float64, 1}, horse_probs::Array{Float64, 1}, p_bird_to_mosquito::Float64, p_mosquito_to_bird::Float64, p_mosquito_to_human::Float64, p_mosquito_to_horse::Float64; seed_birds::Int=1, seed_mosquitoes::Int=0, seed_humans::Int=0, seed_horses::Int=0, bird_incubation_time::Int=10, mosquito_incubation_time::Int=4, human_prob::Union{Float64, AbstractVector{<:Float64}, Nothing}=nothing)

  status_birds = zeros(Int, n_birds)
  latent_birds = zeros(Int, n_birds)

  status_mosquitoes = zeros(Int, n_mosquitoes)
  latent_mosquitoes = zeros(Int, n_mosquitoes)
  status_humans = zeros(Int, n_humans)
  status_horses = zeros(Int, n_horses)

  infected_ever = falses(n_mosquitoes)
  mosq_ever_total = 0

  age_mosquitoes = Random.rand(Distributions.DiscreteUniform(0, mosquito_life_span), n_mosquitoes)

  human_weight = if human_prob === nothing
    n_humans > 0 ? ones(Float64, n_humans) : Float64[]
  elseif human_prob isa Number
    fill(Float64(human_prob), n_humans)
  else
    Float64.(human_prob)
  end
  if length(human_weight) != n_humans
    error("human_prob must have length n_humans when provided as a vector")
  end

  host_weight_bird = n_birds > 0 ? 1.0 : 0.0
  host_weight_human = n_humans > 0 ? Statistics.mean(human_weight) : 0.0
  host_weight_horse = n_horses > 0 ? 1.0 : 0.0

  if human_prob === nothing
    bird_group_weight = 1.0
    human_group_weight = 1.0
    horse_group_weight = 1.0
  else
    total_host_weight = host_weight_bird + host_weight_human + host_weight_horse
    bird_group_weight = total_host_weight > 0 ? host_weight_bird / total_host_weight : 0.0
    human_group_weight = total_host_weight > 0 ? host_weight_human / total_host_weight : 0.0
    horse_group_weight = total_host_weight > 0 ? host_weight_horse / total_host_weight : 0.0
  end

  if seed_birds > 0 && n_birds > 0
    status_birds[Random.rand(1:n_birds, min(seed_birds, n_birds))] .= 1
  end
  if seed_mosquitoes > 0 && n_mosquitoes > 0
    seeded_idx = Random.rand(1:n_mosquitoes, min(seed_mosquitoes, n_mosquitoes))
    status_mosquitoes[seeded_idx] .= 1
    infected_ever[seeded_idx] .= true
    mosq_ever_total += length(seeded_idx)
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

  n_bird_infections[1] = sum(status_birds .> 0) + sum(latent_birds .> 0)
  n_mosquito_infections[1] = sum(status_mosquitoes .> 0) + sum(latent_mosquitoes .> 0)
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
    dead_ms = findall(age_mosquitoes .> mosquito_life_span)
    if !isempty(dead_ms)
      status_mosquitoes[dead_ms] .= 0
      latent_mosquitoes[dead_ms] .= 0
      age_mosquitoes[dead_ms] .= 0
    end

    expiring_birds = findall(latent_birds .== 1)
    latent_birds[findall(latent_birds .> 0)] .-= 1
    if !isempty(expiring_birds)
      status_birds[expiring_birds] .= 1
    end

    expiring_ms = findall(latent_mosquitoes .== 1)
    latent_mosquitoes[findall(latent_mosquitoes .> 0)] .-= 1
    if !isempty(expiring_ms)
      status_mosquitoes[expiring_ms] .= 1
    end

    i_bs = findall(status_birds .> 0)
    s_ms = findall((status_mosquitoes .== 0) .& (latent_mosquitoes .== 0))
    if !isempty(i_bs) && !isempty(s_ms)
      new_mosquito_infections = one_way_bites(bird_probs[i_bs], mosquito_probs[s_ms], p_bird_to_mosquito)
      if sum(new_mosquito_infections) > 0
        if mosquito_incubation_time > 0
          latent_mosquitoes[s_ms[new_mosquito_infections]] .= mosquito_incubation_time
        else
          status_mosquitoes[s_ms[new_mosquito_infections]] .= 1
        end
      end
    end

    i_ms = findall(status_mosquitoes .> 0)

    s_bs = findall((status_birds .== 0) .& (latent_birds .== 0))
    if !isempty(i_ms) && !isempty(s_bs)
      new_bird_infections = one_way_bites(mosquito_probs[i_ms], bird_probs[s_bs] .* bird_group_weight, p_mosquito_to_bird)
      if sum(new_bird_infections) > 0
        if bird_incubation_time > 0
          latent_birds[s_bs[new_bird_infections]] .= bird_incubation_time
        else
          status_birds[s_bs[new_bird_infections]] .= 1
        end
      end
    end

    s_hs = findall(status_humans .== 0)
    if !isempty(i_ms) && !isempty(s_hs)
      human_targets = human_probs[s_hs] .* human_group_weight .* (human_weight[s_hs])
      new_human_infections = one_way_bites(mosquito_probs[i_ms], human_targets, p_mosquito_to_human)
      if sum(new_human_infections) > 0
        status_humans[s_hs[new_human_infections]] .= 1
      end
    end

    s_horses = findall(status_horses .== 0)
    if !isempty(i_ms) && !isempty(s_horses)
      new_horse_infections = one_way_bites(mosquito_probs[i_ms], horse_probs[s_horses] .* horse_group_weight, p_mosquito_to_horse)
      if sum(new_horse_infections) > 0
        status_horses[s_horses[new_horse_infections]] .= 1
      end
    end

    n_bird_infections[s] = sum(status_birds .> 0) + sum(latent_birds .> 0)
    n_mosquito_infections[s] = sum(status_mosquitoes .> 0) + sum(latent_mosquitoes .> 0)
    n_human_infections[s] = sum(status_humans .> 0)
    n_horse_infections[s] = sum(status_horses .> 0)

    n_bird_recovered[s] = sum(status_birds .< 0)
    n_human_recovered[s] = sum(status_humans .< 0)
    n_horse_recovered[s] = sum(status_horses .< 0)

  end

  return n_mosquito_infections, n_bird_infections, n_human_infections, n_horse_infections, n_bird_recovered, n_human_recovered, n_horse_recovered

end

"""
    BiteEdge

Element type of the optional `network_edges` log of [`bite_steps_quad_decay`](@ref): one record per bite,
`(step, mosquito, target_type, target, success)`, where `success` is `true` if the bite transmitted WNV in
either direction.
"""
const BiteEdge = NamedTuple{(:step, :mosquito, :target_type, :target, :success), Tuple{Int, Int, Symbol, Int, Bool}}

"""
    TransmissionEvent

Element type of the optional `transmission_log` of [`bite_steps_quad_decay`](@ref): one record per infection,
`(step, source_type, source, target_type, target)`, i.e. who infected whom.
"""
const TransmissionEvent = NamedTuple{(:step, :source_type, :source, :target_type, :target), Tuple{Int, Symbol, Int, Symbol, Int}}

export BiteEdge, TransmissionEvent

# Cumulative per-individual weights of one host type, for drawing an individual with probability
# proportional to its weight. Uniform weights take the O(1) path.
struct HostWeights
  n::Int
  total::Float64
  uniform::Bool
  cum::Vector{Float64}
end

function HostWeights(w::AbstractVector{<:Real}, n::Int, name::AbstractString)
  n == 0 && return HostWeights(0, 0.0, true, Float64[])
  length(w) == n || throw(ArgumentError("$(name) has length $(length(w)) but the population size is $(n)"))
  any(x -> x < 0, w) && throw(ArgumentError("$(name) must be non-negative"))
  cum = cumsum(Float64.(w))
  total = cum[end]
  uniform = total > 0 && all(x -> x == w[1], w)
  return HostWeights(n, total, uniform, cum)
end

@inline function pick_individual(hw::HostWeights)
  hw.uniform && return rand(1:hw.n)
  idx = searchsortedlast(hw.cum, rand() * hw.total) + 1
  return min(idx, hw.n)
end

# counts[g + 1] += 1, growing the vector as needed (generation 0 lives at index 1)
@inline function count_generation!(counts::Vector{Int}, g::Int)
  while length(counts) <= g
    push!(counts, 0)
  end
  counts[g + 1] += 1
  return counts
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
    p_extra_bite::Float64=0.0,
    collect_bite_counts::Bool=false,
    bite_cycle_counts::Dict{Int, Int}=Dict{Int, Int}(),
    network_edges::Union{Nothing, Vector{BiteEdge}}=nothing,
    transmission_log::Union{Nothing, Vector{TransmissionEvent}}=nothing,
    bird_incubation_time::Int=10,
    mosquito_incubation_time::Int=4,
    human_prob=nothing,
  )

Quadripartite West Nile virus simulation (birds = reservoir, mosquitoes = vector, humans and horses =
dead-end hosts) with an explicit mosquito gonotrophic cycle and independent host choice for every bite.

## Biting model
Each mosquito has a gonotrophic cycle of `gonotrophic_length` steps with a random phase (so feeding is spread
evenly over steps, and newborn mosquitoes draw a fresh phase). On the first step of each cycle the mosquito takes
one blood meal. A meal consists of one bite plus, after every bite, a further bite with probability
`p_extra_bite` (so bites per meal are geometric with mean `1 / (1 - p_extra_bite)`; `p_extra_bite` is the
quantity that mixed blood-meal data estimate). Every bite chooses its host independently of all other bites and
mosquitoes: host type with probability proportional to the summed weights of that type (`sum(bird_probs)`,
`sum(human_probs)`, `sum(horse_probs)`), and the individual within the type with probability proportional to
its own weight. With uniform weights the share of bites that land on birds is therefore
`n_birds * bird_probs[1] / (n_birds * bird_probs[1] + n_humans * human_probs[1] + n_horses * horse_probs[1])`.
`mosquito_probs` is accepted for signature compatibility but not used (all mosquitoes bite at the same rate).

## Infection model
Newly infected birds and mosquitoes are latent for `bird_incubation_time` / `mosquito_incubation_time` steps,
then infectious: birds for `bird_infection_time` steps, mosquitoes until death. Humans and horses are infected
for `*_infection_time` steps and never transmit. Mosquitoes live `mosquito_life_span + 1` steps and are replaced
by susceptible newborns. Seeded birds become infectious on the first biting step (step 2) and get the full
infectious period, so the index case is comparable with later cases.

## Returns
A `NamedTuple`; the first eight fields are, in order, the values returned by earlier versions, so
`mosq, bird, hum, horse, bird_rec, hum_rec, horse_rec, mosq_ever = bite_steps_quad_decay(...)` still works:

- `n_mosquito_infections`, `n_bird_infections`, `n_human_infections`, `n_horse_infections`: infected (latent +
  infectious) counts per step;
- `n_bird_recovered`, `n_human_recovered`, `n_horse_recovered`: recovered counts per step;
- `mosq_ever_total`: number of individual mosquitoes infected over the run (a replacement newborn is a new
  individual);
- `bird_generations`, `mosquito_generations`, `human_generations`, `horse_generations`: number of infections
  by transmission generation, `v[g + 1]` being generation `g`. Seeded individuals are generation 0; a mosquito
  infected by a generation-`g` bird is generation `g`; a host infected by a generation-`g` mosquito is
  generation `g + 1`. With one seeded bird, `bird_generations[2]` is the realised host-to-host reproduction
  number of the index case and `mosquito_generations[1]` the number of mosquitoes it infected; averaged over
  replicates these estimate `R0` and `R_hv` of [`r0_analytic`](@ref);
- `mosq_individuals_total`: mosquito individuals that lived during the run (initial population plus
  replacements), the denominator for a per-individual mosquito attack rate.

## Optional logging
- `collect_bite_counts=true` with a `Dict{Int,Int}` in `bite_cycle_counts`: histogram of bites per blood meal;
- `network_edges`: a `Vector{BiteEdge}` that receives one record per bite;
- `transmission_log`: a `Vector{TransmissionEvent}` that receives one record per infection (who infected whom).
All three containers are emptied at the start of the run.

## Migration from versions before 0.2
`bite_decay` has been removed. Earlier versions walked every mosquito through one shared random host order per
step, which concentrated bites on the first hosts in that order, and the chance of extra bites grew with the
number of hosts. Replace `bite_decay=δ` with `p_extra_bite=q`, where `q` is the probability of a further bite
after each bite (`q = 0` means exactly one bite per cycle).
"""
function bite_steps_quad_decay(n_steps::Int64, n_birds::Int64, n_mosquitoes::Int64, n_humans::Int64, n_horses::Int64, bird_infection_time::Int64, human_infection_time::Int64, horse_infection_time::Int64, mosquito_life_span::Int64, bird_probs::Array{Float64, 1}, mosquito_probs::Array{Float64, 1}, human_probs::Array{Float64, 1}, horse_probs::Array{Float64, 1}, p_bird_to_mosquito::Float64, p_mosquito_to_bird::Float64, p_mosquito_to_human::Float64, p_mosquito_to_horse::Float64; seed_birds::Int=1, seed_mosquitoes::Int=0, seed_humans::Int=0, seed_horses::Int=0, gonotrophic_length::Int=4, p_extra_bite::Float64=0.0, collect_bite_counts::Bool=false, bite_cycle_counts::Dict{Int, Int}=Dict{Int, Int}(), network_edges::Union{Nothing, Vector{BiteEdge}}=nothing, transmission_log::Union{Nothing, Vector{TransmissionEvent}}=nothing, bird_incubation_time::Int=10, mosquito_incubation_time::Int=4, human_prob::Union{Float64, AbstractVector{<:Float64}, Nothing}=nothing, bite_decay=nothing)

  bite_decay === nothing || throw(ArgumentError("`bite_decay` has been removed; use `p_extra_bite` (probability of a further bite after each bite within a gonotrophic cycle). See the `bite_steps_quad_decay` docstring."))
  n_steps >= 1 || throw(ArgumentError("n_steps must be at least 1"))
  gonotrophic_length >= 1 || throw(ArgumentError("gonotrophic_length must be at least 1"))
  (0.0 <= p_extra_bite < 1.0) || throw(ArgumentError("p_extra_bite must be in [0, 1)"))

  birds_w = HostWeights(bird_probs, n_birds, "bird_probs")
  humans_w = HostWeights(human_probs, n_humans, "human_probs")
  horses_w = HostWeights(horse_probs, n_horses, "horse_probs")
  total_w = birds_w.total + humans_w.total + horses_w.total
  bird_cut = birds_w.total
  human_cut = birds_w.total + humans_w.total

  # status: 0 = susceptible, >0 = infectious (steps since becoming infectious), <0 = recovered
  status_birds = zeros(Int, n_birds)
  latent_birds = zeros(Int, n_birds)
  status_mosquitoes = zeros(Int, n_mosquitoes)
  latent_mosquitoes = zeros(Int, n_mosquitoes)
  status_humans = zeros(Int, n_humans)
  status_horses = zeros(Int, n_horses)

  age_mosquitoes = Random.rand(Distributions.DiscreteUniform(0, mosquito_life_span), n_mosquitoes)
  L = gonotrophic_length
  cycle_day = Random.rand(0:(L - 1), n_mosquitoes)   # phase within the gonotrophic cycle; feeds when it wraps to 0

  # transmission generation of each infected individual (-1 = never infected)
  gen_birds = fill(-1, n_birds)
  gen_mosquitoes = fill(-1, n_mosquitoes)
  bird_generations = Int[]
  mosquito_generations = Int[]
  human_generations = Int[]
  horse_generations = Int[]
  mosq_ever_total = 0
  mosq_births = 0

  if seed_birds > 0 && n_birds > 0
    for idx in Random.randperm(n_birds)[1:min(seed_birds, n_birds)]
      latent_birds[idx] = 1      # infectious from step 2 with the full infectious period
      gen_birds[idx] = 0
      count_generation!(bird_generations, 0)
    end
  end
  if seed_mosquitoes > 0 && n_mosquitoes > 0
    for m in Random.randperm(n_mosquitoes)[1:min(seed_mosquitoes, n_mosquitoes)]
      status_mosquitoes[m] = 1
      gen_mosquitoes[m] = 0
      count_generation!(mosquito_generations, 0)
      mosq_ever_total += 1
    end
  end
  if seed_humans > 0 && n_humans > 0
    for idx in Random.randperm(n_humans)[1:min(seed_humans, n_humans)]
      status_humans[idx] = 1
      count_generation!(human_generations, 0)
    end
  end
  if seed_horses > 0 && n_horses > 0
    for idx in Random.randperm(n_horses)[1:min(seed_horses, n_horses)]
      status_horses[idx] = 1
      count_generation!(horse_generations, 0)
    end
  end

  collect_bite_counts && empty!(bite_cycle_counts)
  network_edges === nothing || empty!(network_edges)
  transmission_log === nothing || empty!(transmission_log)

  n_bird_infections = Vector{Int64}(undef, n_steps)
  n_mosquito_infections = Vector{Int64}(undef, n_steps)
  n_human_infections = Vector{Int64}(undef, n_steps)
  n_horse_infections = Vector{Int64}(undef, n_steps)
  n_bird_recovered = Vector{Int64}(undef, n_steps)
  n_human_recovered = Vector{Int64}(undef, n_steps)
  n_horse_recovered = Vector{Int64}(undef, n_steps)

  n_bird_infections[1] = count(>(0), status_birds) + count(>(0), latent_birds)
  n_mosquito_infections[1] = count(>(0), status_mosquitoes) + count(>(0), latent_mosquitoes)
  n_human_infections[1] = count(>(0), status_humans)
  n_horse_infections[1] = count(>(0), status_horses)
  n_bird_recovered[1] = count(<(0), status_birds)
  n_human_recovered[1] = count(<(0), status_humans)
  n_horse_recovered[1] = count(<(0), status_horses)

  for s = 2:n_steps

    # infectious periods run their course
    @inbounds for i in 1:n_birds
      if status_birds[i] > 0
        status_birds[i] += 1
        status_birds[i] > bird_infection_time && (status_birds[i] = -1)
      end
    end
    @inbounds for i in 1:n_humans
      if status_humans[i] > 0
        status_humans[i] += 1
        status_humans[i] > human_infection_time && (status_humans[i] = -1)
      end
    end
    @inbounds for i in 1:n_horses
      if status_horses[i] > 0
        status_horses[i] += 1
        status_horses[i] > horse_infection_time && (status_horses[i] = -1)
      end
    end

    # mosquitoes age; those past the life span are replaced by susceptible newborns with a fresh cycle phase
    @inbounds for m in 1:n_mosquitoes
      age_mosquitoes[m] += 1
      if age_mosquitoes[m] > mosquito_life_span
        age_mosquitoes[m] = 0
        mosq_births += 1
        status_mosquitoes[m] = 0
        latent_mosquitoes[m] = 0
        gen_mosquitoes[m] = -1
        cycle_day[m] = rand(0:(L - 1))
      end
    end

    # latent individuals become infectious
    @inbounds for i in 1:n_birds
      if latent_birds[i] == 1
        latent_birds[i] = 0
        status_birds[i] = 1
      elseif latent_birds[i] > 1
        latent_birds[i] -= 1
      end
    end
    @inbounds for m in 1:n_mosquitoes
      if latent_mosquitoes[m] == 1
        latent_mosquitoes[m] = 0
        status_mosquitoes[m] = 1
      elseif latent_mosquitoes[m] > 1
        latent_mosquitoes[m] -= 1
      end
    end

    # blood meals: a mosquito whose cycle wraps around this step feeds; each bite picks its own host
    @inbounds for m in 1:n_mosquitoes
      cycle_day[m] += 1
      cycle_day[m] < L && continue
      cycle_day[m] = 0
      total_w > 0 || continue
      n_bites = 0
      while true
        n_bites += 1
        u = rand() * total_w
        success = false
        if u < bird_cut
          idx = pick_individual(birds_w)
          if status_mosquitoes[m] > 0 && status_birds[idx] == 0 && latent_birds[idx] == 0
            success = rand() < p_mosquito_to_bird
            if success
              if bird_incubation_time > 0
                latent_birds[idx] = bird_incubation_time
              else
                status_birds[idx] = 1
              end
              g = gen_mosquitoes[m] + 1
              gen_birds[idx] = g
              count_generation!(bird_generations, g)
              transmission_log === nothing || push!(transmission_log, (step=s, source_type=:mosquito, source=m, target_type=:bird, target=idx))
            end
          elseif status_birds[idx] > 0 && status_mosquitoes[m] == 0 && latent_mosquitoes[m] == 0
            success = rand() < p_bird_to_mosquito
            if success
              if mosquito_incubation_time > 0
                latent_mosquitoes[m] = mosquito_incubation_time
              else
                status_mosquitoes[m] = 1
              end
              mosq_ever_total += 1
              g = gen_birds[idx]
              gen_mosquitoes[m] = g
              count_generation!(mosquito_generations, g)
              transmission_log === nothing || push!(transmission_log, (step=s, source_type=:bird, source=idx, target_type=:mosquito, target=m))
            end
          end
          network_edges === nothing || push!(network_edges, (step=s, mosquito=m, target_type=:bird, target=idx, success=success))
        elseif u < human_cut
          idx = pick_individual(humans_w)
          if status_mosquitoes[m] > 0 && status_humans[idx] == 0
            success = rand() < p_mosquito_to_human
            if success
              status_humans[idx] = 1
              count_generation!(human_generations, gen_mosquitoes[m] + 1)
              transmission_log === nothing || push!(transmission_log, (step=s, source_type=:mosquito, source=m, target_type=:human, target=idx))
            end
          end
          network_edges === nothing || push!(network_edges, (step=s, mosquito=m, target_type=:human, target=idx, success=success))
        else
          idx = pick_individual(horses_w)
          if status_mosquitoes[m] > 0 && status_horses[idx] == 0
            success = rand() < p_mosquito_to_horse
            if success
              status_horses[idx] = 1
              count_generation!(horse_generations, gen_mosquitoes[m] + 1)
              transmission_log === nothing || push!(transmission_log, (step=s, source_type=:mosquito, source=m, target_type=:horse, target=idx))
            end
          end
          network_edges === nothing || push!(network_edges, (step=s, mosquito=m, target_type=:horse, target=idx, success=success))
        end
        rand() < p_extra_bite || break
      end
      if collect_bite_counts
        bite_cycle_counts[n_bites] = get(bite_cycle_counts, n_bites, 0) + 1
      end
    end

    n_bird_infections[s] = count(>(0), status_birds) + count(>(0), latent_birds)
    n_mosquito_infections[s] = count(>(0), status_mosquitoes) + count(>(0), latent_mosquitoes)
    n_human_infections[s] = count(>(0), status_humans)
    n_horse_infections[s] = count(>(0), status_horses)
    n_bird_recovered[s] = count(<(0), status_birds)
    n_human_recovered[s] = count(<(0), status_humans)
    n_horse_recovered[s] = count(<(0), status_horses)

  end

  return (
    n_mosquito_infections = n_mosquito_infections,
    n_bird_infections = n_bird_infections,
    n_human_infections = n_human_infections,
    n_horse_infections = n_horse_infections,
    n_bird_recovered = n_bird_recovered,
    n_human_recovered = n_human_recovered,
    n_horse_recovered = n_horse_recovered,
    mosq_ever_total = mosq_ever_total,
    bird_generations = bird_generations,
    mosquito_generations = mosquito_generations,
    human_generations = human_generations,
    horse_generations = horse_generations,
    mosq_individuals_total = n_mosquitoes + mosq_births,
  )

end

export r0_analytic
"""
    r0_analytic(; n_birds, n_mosquitoes, n_humans=0, n_horses=0,
                bird_weight, human_weight=0.0, horse_weight=0.0,
                p_bird_to_mosquito, p_mosquito_to_bird, p_mosquito_to_human=0.0, p_mosquito_to_horse=0.0,
                bird_infection_time, mosquito_life_span, mosquito_incubation_time,
                gonotrophic_length, p_extra_bite=0.0)

Expected next-generation quantities of [`bite_steps_quad_decay`](@ref) in a fully susceptible population
(Ross–Macdonald-type calculation that follows the discrete-step mechanics of the simulation). Weights may be
scalars (uniform within the type) or the per-individual vectors passed to the simulation.

Returns a `NamedTuple` with
- `R_hv`: mosquitoes infected by one infectious bird;
- `R_vh`: birds infected by one infected mosquito (averaged over its age at infection);
- `R0 = R_hv * R_vh`: host-to-host basic reproduction number, which `bird_generations[2]` of a single-seed
  simulation estimates; `R0_geometric = sqrt(R0)` is the per-generation (host→vector→host) convention;
- `spillover_human`, `spillover_horse`: humans/horses infected by the mosquitoes that one infectious bird infects;
- `share_bird`, `share_human`, `share_horse`: proportion of bites landing on each host type;
- `bites_per_mosquito_per_step`, `infectious_meals_per_infected_mosquito`.

Assumes `mosquito_incubation_time >= 1` (with no incubation a mosquito can transmit within the meal that
infected it, which this calculation does not count).
"""
function r0_analytic(; n_birds::Int, n_mosquitoes::Int, n_humans::Int=0, n_horses::Int=0,
    bird_weight, human_weight=0.0, horse_weight=0.0,
    p_bird_to_mosquito::Real, p_mosquito_to_bird::Real, p_mosquito_to_human::Real=0.0, p_mosquito_to_horse::Real=0.0,
    bird_infection_time::Int, mosquito_life_span::Int, mosquito_incubation_time::Int,
    gonotrophic_length::Int, p_extra_bite::Real=0.0)

  total_weight(w::Real, n::Int) = n > 0 ? Float64(w) * n : 0.0
  total_weight(w::AbstractVector, n::Int) = n > 0 ? Float64(sum(w)) : 0.0

  W_b = total_weight(bird_weight, n_birds)
  W_h = total_weight(human_weight, n_humans)
  W_o = total_weight(horse_weight, n_horses)
  W = W_b + W_h + W_o
  share_b = W > 0 ? W_b / W : 0.0
  share_h = W > 0 ? W_h / W : 0.0
  share_o = W > 0 ? W_o / W : 0.0

  L = max(gonotrophic_length, 1)
  bites_per_meal = 1 / (1 - p_extra_bite)
  bites_per_step = bites_per_meal / L

  # bites received by one bird per step, and mosquitoes infected by one bird over its infectious period
  bites_on_bird = n_birds > 0 ? bites_per_step * share_b * n_mosquitoes / n_birds : 0.0
  R_hv = bird_infection_time * bites_on_bird * p_bird_to_mosquito

  # a mosquito is infected on a feeding day at age a ~ Uniform(0..life_span); it feeds again at ages a + kL,
  # is infectious once k*L >= incubation, and bites while its age is <= life_span
  Λ = mosquito_life_span
  k_min = cld(mosquito_incubation_time, L)
  meals = 0.0
  for a in 0:Λ
    meals += max(0, fld(Λ - a, L) - k_min + 1)
  end
  meals /= (Λ + 1)

  R_vh = meals * bites_per_meal * share_b * p_mosquito_to_bird
  R0 = R_hv * R_vh
  spill_h = R_hv * meals * bites_per_meal * share_h * p_mosquito_to_human
  spill_o = R_hv * meals * bites_per_meal * share_o * p_mosquito_to_horse

  return (
    R_hv = R_hv,
    R_vh = R_vh,
    R0 = R0,
    R0_geometric = sqrt(R0),
    spillover_human = spill_h,
    spillover_horse = spill_o,
    share_bird = share_b,
    share_human = share_h,
    share_horse = share_o,
    bites_per_mosquito_per_step = bites_per_step,
    infectious_meals_per_infected_mosquito = meals,
  )
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

!!! warning "Deprecated"
    This function does not compute a reproduction number. With `seed_cases > 0` it returns, per replicate, the
    smallest infected count that ever exceeds the seed (typically 2 or 3), and with `seed_cases == 0` the ratio
    `(second-smallest positive count - smallest positive count) / smallest positive count`. It is kept for
    backward compatibility and will be removed. Use the `bird_generations` / `mosquito_generations` outputs of
    [`bite_steps_quad_decay`](@ref) (realised reproduction numbers) or [`r0_analytic`](@ref) instead.

Returns a named Tuple containing (1) the mean over replicates of the quantity described above, (2) its value for each replicate, (3) its running mean over successive numbers of replicates.

## Parameters
* `n_reps` The number of repetitions in the simulation. Int64. 
* `infection_ts` The matrix time series of number of infected cases at each time step. A two-dimensional array of Int64 with dimensions equal to number of repetitions and number of time steps.
* `seed_cases` Number of infected cases seeded into the population at the beginning of each repetition. If `seed_cases=0`, then the first cases will be detected in the data. (For example, if we seed 1 infected human case into the population, then we would use `seed_cases=1` for the human R0 and `seed_cases=0` for the mosquito R0.)
"""
function calculate_r0(n_reps::Int64, infection_ts::Array{Int64,2}, seed_cases::Int64)

  Base.depwarn("calculate_r0 does not compute a reproduction number; use the generation counts returned by bite_steps_quad_decay or r0_analytic", :calculate_r0)

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
