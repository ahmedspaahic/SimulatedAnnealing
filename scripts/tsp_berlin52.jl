# Simulirano kaljenje na Berlin52 instanci problema trgovačkog putnika (TSPLIB).
# Poznata optimalna dužina ture: 7542.
#
# Pokretanje (iz korijena repozitorija):
#     julia --project=. scripts/tsp_berlin52.jl

include(joinpath(@__DIR__, "..", "src", "SimulatedAnnealing.jl"))
using .SimulatedAnnealing

using DelimitedFiles
using Plots
using Random

Random.seed!(42)  # za ponovljive rezultate

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")
const OPTIMUM     = 7542.0
mkpath(RESULTS_DIR)

# =========================
# Podaci za Berlin52
# =========================
cities = readdlm(joinpath(@__DIR__, "..", "data", "berlin52.txt"), Float64; skipstart=1)

# =========================
# Pokretanje SA 30 puta i čuvanje najbolje ture
# =========================
function main(cities; num_runs=30)
    best_tour = Int[]
    best_len  = Inf

    for run in 1:num_runs
        tour, len = tsp_sa(cities)
        println("Run $run: dužina ture = ", round(len, digits=2))
        if len < best_len
            best_len  = len
            best_tour = tour
        end
    end

    rel_err = (best_len - OPTIMUM) / OPTIMUM * 100
    println("\nNajbolja dužina ture nakon $num_runs pokretanja = ", round(best_len, digits=2))
    println("Relativna greška u odnosu na optimum ($OPTIMUM) = ", round(rel_err, digits=2), " %")

    # =========================
    # Crtanje najbolje ture
    # =========================
    x = cities[:,1]
    y = cities[:,2]
    tour = [best_tour; best_tour[1]]  # zatvori turu

    p = plot(x[tour], y[tour], seriestype=:path, marker=:circle,
             legend=false,
             title="Najbolja dužina ture = $(round(best_len, digits=2))",
             xlabel="X", ylabel="Y")
    savefig(p, joinpath(RESULTS_DIR, "tsp_berlin52.png"))
    println("Plot sačuvan u ", normpath(RESULTS_DIR))
end

main(cities)
