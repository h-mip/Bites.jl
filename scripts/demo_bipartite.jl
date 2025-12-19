using Bites, Random, Plots

Random.seed!(123)

n_steps = 50
n_humans = 12
n_mosquitoes = 60
human_infection_time = 4
mosquito_life_span = 20
transmission_prob = 0.35

# Uniform bite propensities for simplicity
human_probs = fill(1 / n_mosquitoes, n_humans)
mosquito_probs = fill(1 / n_humans, n_mosquitoes)

mosq_ts, hum_ts, hum_rec = bite_steps(
  n_steps,
  n_humans,
  n_mosquitoes,
  human_infection_time,
  mosquito_life_span,
  human_probs,
  mosquito_probs,
  transmission_prob,
)

println("Mosquito infected counts: ", mosq_ts)
println("Human infected counts:    ", hum_ts)
println("Human recovered counts:   ", hum_rec)

ts = 1:n_steps

p1 = plot(ts, mosq_ts; label="Mosquito infected", lw=2, xlabel="Step", ylabel="Count", title="Mosquito infections")
p2 = plot(ts, hum_ts; label="Human infected", lw=2, xlabel="Step", ylabel="Count", title="Human infections")
p = plot(p1, p2; layout=(2, 1), size=(800, 600), title="Bipartite infections", link=:none)

out_dir = joinpath(@__DIR__, "..", "plots")
mkpath(out_dir)
out_path = joinpath(out_dir, "demo_bipartite.png")
savefig(p, out_path)
println("Saved plot to " * out_path)
