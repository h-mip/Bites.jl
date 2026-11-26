```@meta
CurrentModule = Bites
```

# Bites

[![DOI](https://zenodo.org/badge/398859191.svg)](https://doi.org/10.5281/zenodo.18088794)

Bites is a Julia package for analyzing mosquito-host biting networks. In addition to bipartite mosquito-host models, Bites now supports quadripartite simulations with birds, mosquitoes, humans, and horses, suitable for West Nile Virus (WNV) transmission modeling. The package includes features for simulating incubation periods, customizing mosquito->human targeting weights, and an optional bite-decay variant that reduces a mosquito’s chance of additional bites within a gonotrophic cycle across all host types.

Bites.jl was developed as part of the Host-Mosquito Interaction Project (H-MIP), funded by the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (Grant agreement No. 853271). The quadripartite network components used for West Nile Virus modeling were added as part of the E4Warning project (Eco-Epidemiological Intelligence for early Warning and response to mosquito-borne disease risk in Endemic and Emergence settings) funded by the European Union's Horizon Europe programme under Grant Agreement 101086640.

## Quick start

```julia
using Bites, Random

# simple bipartite example
n_steps = 50
n_humans = 12
n_mosquitoes = 60
human_probs = fill(1 / n_mosquitoes, n_humans)
mosquito_probs = fill(1 / n_humans, n_mosquitoes)
mosq_ts, hum_ts, hum_rec = bite_steps(n_steps, n_humans, n_mosquitoes, 4, 20, human_probs, mosquito_probs, 0.35)

# quadripartite example
bird_probs = fill(0.8, 5_000)
horse_probs = fill(0.6, 300)
mosq_ts, bird_ts, hum_ts, horse_ts, bird_rec, hum_rec, horse_rec = bite_steps_quad(
	50, 5_000, 5_000, 100, 300, 4, 5, 5, 20,
	bird_probs, mosquito_probs, human_probs, horse_probs,
	0.4, 0.6, 0.25, 0.2;
	seed_birds=1,
)
```

## Incubation and human targeting

- Birds and mosquitoes now pass through a latent/incubation period before becoming infectious. Control durations with `bird_incubation_time` and `mosquito_incubation_time` (defaults 10 and 4). Latent individuals are counted in the infection timeseries but do not transmit until their timers expire.
- Human bite weight can be customized via `human_prob` (scalar or vector). This weight is normalized against bird and horse weights so you can emphasize/de-emphasize humans without changing other host probabilities. Leaving it `nothing` preserves the previous behavior.

Example with explicit incubation and human weighting:

```julia
mosq_ts, bird_ts, hum_ts, horse_ts, bird_rec, hum_rec, horse_rec = bite_steps_quad(
	60, 5_000, 5_000, 1_000, 3_000, 4, 5, 5, 20,
	bird_probs, mosquito_probs, human_probs, horse_probs,
	0.3, 0.4, 0.25, 0.2;
	bird_incubation_time=8,
	mosquito_incubation_time=4,
	human_prob=0.5,
	seed_birds=1,
)
```

## Quadripartite model with a gonotrophic cycle: `bite_steps_quad_decay`

`bite_steps_quad_decay` is the recommended WNV model. Each mosquito has a gonotrophic cycle of `gonotrophic_length` steps with a random phase and takes one blood meal per cycle. After every bite it takes a further bite with probability `p_extra_bite`, so bites per meal are geometric with mean `1 / (1 - p_extra_bite)`; this is the quantity that mixed blood-meal data estimate. Every bite chooses its host independently: the host type with probability proportional to the summed weights of that type, and the individual with probability proportional to its own weight, so per-individual heterogeneity in attractiveness can be modelled through the `*_probs` vectors.

Key parameters:
- `gonotrophic_length` (Int, default 4): steps between blood meals.
- `p_extra_bite` (Float64, default 0.0): probability of a further bite after each bite within a meal.

The function returns a `NamedTuple` whose first eight fields are the infection and recovery time series and the number of mosquitoes ever infected, followed by infection counts by transmission generation (`bird_generations`, `mosquito_generations`, `human_generations`, `horse_generations`). With one seeded bird, `bird_generations[2]` is the realised host-to-host reproduction number of the index case. Optional logs: `bite_cycle_counts` (bites per meal), `network_edges` (every bite, as `BiteEdge`), `transmission_log` (every infection, as `TransmissionEvent`).

Versions before 0.2 used a `bite_decay` parameter and one shared host order per step; see `CHANGELOG.md` for why this was replaced.

Example:

```julia
out = bite_steps_quad_decay(
	200, 5_000, 5_000, 1_000, 3_000, 4, 5, 5, 20,
	bird_probs, mosquito_probs, human_probs, horse_probs,
	0.3, 0.4, 0.25, 0.2;
	seed_birds=1,
	gonotrophic_length=4,
	p_extra_bite=0.5,
)
out.n_human_infections    # time series of infected humans
out.bird_generations[2]   # birds infected via the mosquitoes the index bird infected
```

## Probability generation and R₀

- `distribute_bite_probabilities(D1, D2, N1, N2, expected_bites)` to draw human/mosquito bite propensities from distributions in `Distributions.jl` and scale to a target expected bites.
- `r0_analytic(; ...)` gives the expected `R_hv`, `R_vh`, `R0` and spillover per index bird for `bite_steps_quad_decay` in a fully susceptible population, following the simulation's discrete-step mechanics.
- `bite_steps_tracked` runs the bipartite model of `bite_steps` and also returns infection counts by transmission generation, so `human_generations[2]` is the realised reproduction number of the index human. Its `seed` keyword chooses the index human uniformly at random (`:uniform`) or in proportion to bite probability (`:degree`).
- `r0_analytic_bipartite(human_probs, mosquito_probs; ...)` gives the expected reproduction number of the bipartite model in a fully susceptible population for both kinds of index case.
- `calculate_r0` is deprecated: it does not compute a reproduction number (see its docstring).

## API reference

```@docs
bite_steps
bite_steps_tracked
bite_steps_quad
bite_steps_quad_decay
r0_analytic
r0_analytic_bipartite
distribute_bite_probabilities
calculate_r0
```

```@index
```

```@autodocs
Modules = [Bites]
```


## About

Bites.jl was developed as part of the [Host-Mosquito Interaction Project (H-MIP)](https://h-mip.github.io/), funded by the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 853271). The quadripartite network components used for West Nile Virus modeling were added as part of the [E4Warning project](https://www.e4warning.eu) (Eco-Epidemiological Intelligence for early Warning and response to mosquito-borne disease risk in Endemic and Emergence settings) funded by the European Union’s Horizon Europe programme under Grant Agreement 101086640.

Copyright 2021-2025 John R.B. Palmer

Bites.jl is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Bites.jl is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
