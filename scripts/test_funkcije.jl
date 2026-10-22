# Poređenje simuliranog kaljenja (SA) i Nelder-Mead metode (NM) na
# četiri benchmark funkcije, uz contour plotove pronađenih rješenja.
#
# Pokretanje (iz korijena repozitorija):
#     julia --project=. scripts/test_funkcije.jl

include(joinpath(@__DIR__, "..", "src", "SimulatedAnnealing.jl"))
using .SimulatedAnnealing

using DataFrames
using Optim
using PrettyTables
using Plots
using Random

Random.seed!(42)  # za ponovljive rezultate

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")
mkpath(RESULTS_DIR)

function sa_vs_nm_table(f, naziv, opseg, global_val;
                        num_sim=3, T0_list=[1.0,3.0,5.0], alpha_list=[0.995,0.998])

    df = DataFrame(Funkcija=String[], T0=Float64[], Alpha=Float64[],
                   SA_X1=Float64[], SA_X2=Float64[], SA_F=Float64[], SA_Diff=Float64[],
                   NM_X1=Float64[], NM_X2=Float64[], NM_F=Float64[], NM_Diff=Float64[])

    # =========================
    # Random start za Nelder-Mead
    # =========================
    x0_nm = [rand()*(opseg[1][2]-opseg[1][1]) + opseg[1][1],
             rand()*(opseg[2][2]-opseg[2][1]) + opseg[2][1]]

    res_nm = optimize(x -> f(x[1], x[2]),
                      [opseg[1][1], opseg[2][1]],   # donja granica
                      [opseg[1][2], opseg[2][2]],   # gornja granica
                      x0_nm,                        # random start
                      Fminbox(NelderMead()))
    x_nm = Optim.minimizer(res_nm)
    f_nm = Optim.minimum(res_nm)
    diff_nm = abs(f_nm - global_val)

    # =========================
    # Višestruke simulacije SA
    # =========================
    for T0 in T0_list
        for alpha in alpha_list
            for i in 1:num_sim
                x0      = [opseg[1][2]-1, opseg[2][1]+1]  # početna tačka SA
                delta_x = [0.5,0.5]
                N       = 50_000
                eps     = 1e-8

                x_sa, f_sa = SA(f, opseg, x0, delta_x, N, eps, T0; alpha=alpha)
                diff_sa = abs(f_sa - global_val)

                push!(df, (naziv, T0, alpha,
                           x_sa[1], x_sa[2], f_sa, diff_sa,
                           x_nm[1], x_nm[2], f_nm, diff_nm))
            end
        end
    end

    return df
end

# =========================
# Testne funkcije: (funkcija, naziv, opseg, globalni minimumi, vrijednost minimuma)
# =========================
testovi = [
    (paraboloid, "Paraboloid",   [(-5,5), (-5,5)],             [(0,0)],   0.0),
    (rastrigin,  "Rastrigin",    [(-5.12,5.12), (-5.12,5.12)], [(0,0)],   0.0),
    (dropwave,   "Drop-Wave",    [(-5.12,5.12), (-5.12,5.12)], [(0,0)],  -1.0),
    (holder,     "Holder Table", [(-10,10), (-10,10)],
        [(8.055, 9.6646), (-8.055, 9.6646), (8.055, -9.6646), (-8.055, -9.6646)], -19.2085),
]

df_total = vcat([sa_vs_nm_table(f, naziv, opseg, gval)
                 for (f, naziv, opseg, _, gval) in testovi]...)

show(df_total, allcols=true, allrows=true)
println()

# LaTeX tabela (kao u radu)
open(joinpath(RESULTS_DIR, "rezultati.tex"), "w") do io
    pretty_table(io, df_total, backend = :latex)
end

# =========================
# Funkcija za crtanje kontura
# =========================
function contour_sa(f, naziv, opseg, global_mins, sa_sol)
    xs = range(opseg[1][1], opseg[1][2], length=200)
    ys = range(opseg[2][1], opseg[2][2], length=200)
    Z = [f(x,y) for y in ys, x in xs]

    c = contour(xs, ys, Z, c=:grays, linewidth=1, grid=true, title=naziv)
    scatter!(c, first.(global_mins), last.(global_mins),
             markercolor=:red, markershape=:circle, label="Globalni min")
    scatter!(c, [sa_sol[1]], [sa_sol[2]], markercolor=:blue, markershape=:x,
             markersize=8, label="SA rješenje")
    xlabel!(c, "x")
    ylabel!(c, "y")
    return c
end

# =========================
# Multi-plot 2x2 sa najboljim SA rješenjem za svaku funkciju
# =========================
gr()
plotovi = map(testovi) do (f, naziv, opseg, global_mins, _)
    redovi = df_total[df_total.Funkcija .== naziv, :]
    best   = redovi[argmin(redovi.SA_F), :]
    contour_sa(f, naziv, opseg, global_mins, (best.SA_X1, best.SA_X2))
end

p = plot(plotovi..., layout=(2,2), size=(900,800))
savefig(p, joinpath(RESULTS_DIR, "testne_funkcije.png"))
println("Rezultati sačuvani u ", normpath(RESULTS_DIR))
