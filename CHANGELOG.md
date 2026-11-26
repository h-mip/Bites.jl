# Changelog

## Unreleased (breaking)

### Fixed
- `bite_steps_quad_decay` chose hosts by walking every mosquito through one shared random host order per
  step. The first hosts in that order therefore received a large share of all bites taken in the step
  (hundreds of bites on a single host with 1000 mosquitoes), the gonotrophic cycle never limited feeding
  (all mosquitoes fed in phase), and the number of extra bites per cycle grew with the number of hosts.
  Attack rates and outbreak probabilities were roughly halved as a result. Every bite now chooses its host
  independently (host type with probability proportional to summed weights, individual proportional to its
  own weight), mosquitoes feed once per gonotrophic cycle with random phases, and extra bites follow
  `p_extra_bite`.
- Seeded birds now get the full infectious period (previously one biting step less than later cases).
- `mosq_ever_total` counts individual mosquitoes (a replacement newborn is a new individual), not slots.

### Added
- `bite_steps_tracked`: the bipartite model of `bite_steps` with infection counts by transmission generation
  (`human_generations[2]` is the realised reproduction number of the index human) and a `seed` keyword
  (`:uniform` or `:degree`) for how the index human is chosen. `bite_steps` itself is unchanged.
- `r0_analytic_bipartite`: expected reproduction numbers for the bipartite model under uniform and
  degree-proportional index cases.
- `p_extra_bite` keyword: probability of a further bite after each bite within a meal, so bites per meal are
  geometric with mean `1 / (1 - p_extra_bite)`. Replaces `bite_decay`.
- Per-individual host weights are honoured (previously only the first element of each `*_probs` vector was
  used).
- Infection counts by transmission generation (`bird_generations`, `mosquito_generations`,
  `human_generations`, `horse_generations`) in the return value. With one seeded bird, `bird_generations[2]`
  is the realised host-to-host reproduction number of the index case.
- `mosq_individuals_total` in the return value: mosquitoes that lived during the run (initial population
  plus replacements), the denominator for a per-individual mosquito attack rate.
- Optional `transmission_log` (who infected whom) and the exported element types `BiteEdge` and
  `TransmissionEvent`.
- `r0_analytic`: expected `R_hv`, `R_vh`, `R0` and spillover per index bird under the same mechanics, for
  validation and experimental design.
- `bite_steps_quad_decay` returns a `NamedTuple`; the first eight fields keep their previous order, so
  positional destructuring of earlier versions still works.

### Deprecated
- `calculate_r0`: it returns the smallest count above the seed (or a ratio of the two smallest positive
  counts), not a reproduction number. Use the generation counts or `r0_analytic`. It still works and will be
  removed in a later release.

### Removed
- `bite_decay` keyword of `bite_steps_quad_decay` (passing it raises an `ArgumentError` naming the
  replacement).
- `Plots` as a package dependency (unused by `src/`; the demo scripts still use it).
