```@meta
CurrentModule = Bites
```

# Bites

Quadripartite and bipartite vector–host simulations (birds ↔ mosquitoes, mosquitoes → humans/horses) with helpers for probability generation and summary metrics.

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
