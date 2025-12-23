```@meta
CurrentModule = Bites
```

# Bites

Quadripartite and bipartite vector–host simulations (birds ↔ mosquitoes, humans ↔ mosquitoes, mosquitoes → humans/horses) with helpers for probability generation and summary metrics.

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

## New: bite decay variant (non-breaking)

`bite_steps_quad_decay` is an optional quadripartite variant that tapers each infected mosquito’s chance of taking additional bites within a day. Multiple bites are still possible, but each successive bite is less likely. Existing functions remain unchanged.

Key parameters:
- `gonotrophic_length` (Int, default 4): interprets steps as days; shorter length implies more bite attempts per day.
- `bite_decay` (Float64, default 0.2): per-bite multiplier applied after each bite attempt; smaller values reduce the chance of second/third bites.

Example:

```julia
mosq_ts, bird_ts, hum_ts, horse_ts, bird_rec, hum_rec, horse_rec = bite_steps_quad_decay(
	200, 5_000, 5_000, 1_000, 3_000, 4, 5, 5, 20,
	bird_probs, mosquito_probs, human_probs, horse_probs,
	0.3, 0.4, 0.25, 0.2;
	seed_birds=1,
	gonotrophic_length=4,
	bite_decay=0.2,
)
```

## Probability generation and R₀

- `distribute_bite_probabilities(D1, D2, N1, N2, expected_bites)` to draw human/mosquito bite propensities from distributions in `Distributions.jl` and scale to a target expected bites.
- `calculate_r0(n_reps, infection_ts, seed_cases)` computes mean R₀ and per-rep R₀ series from infection time series.

## API reference

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