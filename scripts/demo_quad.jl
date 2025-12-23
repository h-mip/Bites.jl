using Bites, Random, Plots

Random.seed!(123)

n_steps = 50
n_birds = 5000
n_mosquitoes = 5000
n_humans = 100
n_horses = 300

bird_infection_time = 4
human_infection_time = 5
horse_infection_time = 5
mosquito_life_span = 20

# Simple bite propensity choices for illustration
bird_probs = fill(0.8, n_birds)
mosquito_probs = fill(1 / n_humans, n_mosquitoes)
human_probs = fill(1 / n_mosquitoes, n_humans)
horse_probs = fill(0.6, n_horses)

p_bird_to_mosquito = 0.4
p_mosquito_to_bird = 0.6
p_mosquito_to_human = 0.25
p_mosquito_to_horse = 0.2

mosq_ts, bird_ts, hum_ts, horse_ts, bird_rec, hum_rec, horse_rec = bite_steps_quad(
  n_steps,
  n_birds,
  n_mosquitoes,
  n_humans,
  n_horses,
  bird_infection_time,
  human_infection_time,
  horse_infection_time,
  mosquito_life_span,
  bird_probs,
  mosquito_probs,
  human_probs,
  horse_probs,
  p_bird_to_mosquito,
  p_mosquito_to_bird,
  p_mosquito_to_human,
  p_mosquito_to_horse;
  seed_birds=1,
  seed_mosquitoes=0,
  seed_humans=0,
  seed_horses=0,
)

println("Mosquito infected counts: ", mosq_ts)
println("Bird infected counts:     ", bird_ts)
println("Human infected counts:    ", hum_ts)
println("Horse infected counts:    ", horse_ts)
println("Bird recovered counts:    ", bird_rec)
println("Human recovered counts:   ", hum_rec)
println("Horse recovered counts:   ", horse_rec)

ts = 1:n_steps

p1 = plot(ts, mosq_ts; label="Mosquito infected", lw=2, xlabel="Step", ylabel="Count", title="Mosquito infections")
p2 = plot(ts, bird_ts; label="Bird infected", lw=2, xlabel="Step", ylabel="Count", title="Bird infections")
p3 = plot(ts, hum_ts; label="Human infected", lw=2, xlabel="Step", ylabel="Count", title="Human infections")
p4 = plot(ts, horse_ts; label="Horse infected", lw=2, xlabel="Step", ylabel="Count", title="Horse infections")
p = plot(p1, p2, p3, p4; layout=(2, 2), size=(950, 750), title="Quadripartite infections", link=:none)

out_dir = joinpath(@__DIR__, "..", "plots")
mkpath(out_dir)
out_path = joinpath(out_dir, "demo_quad.png")
savefig(p, out_path)
println("Saved plot to " * out_path)
