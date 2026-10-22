module SimulatedAnnealing

using LinearAlgebra
using Random

export SA, paraboloid, rastrigin, dropwave, holder,
       tour_length, tsp_sa

# =========================
# Simulated Annealing
# =========================
"""
    SA(f, opseg, x0, delta_x, N, eps, T0; alpha=0.995)

Traži globalni minimum funkcije dvije promjenljive `f(x1, x2)` metodom
simuliranog kaljenja.

- `f`       -- funkcija cilja oblika `f(x1, x2)`
- `opseg`   -- opseg pretraživanja oblika `[(xmin, xmax), (ymin, ymax)]`
- `x0`      -- početna tačka pretraživanja
- `delta_x` -- maksimalni pomak po svakoj koordinati
- `N`       -- maksimalan broj iteracija
- `eps`     -- minimalna promjena vrijednosti funkcije između dvije iteracije
- `T0`      -- početna temperatura
- `alpha`   -- faktor geometrijskog hlađenja, `T_{k+1} = alpha * T_k`

Vraća par `(x_best, f_best)`.
"""
function SA(f, opseg, x0, delta_x, N, eps, T0; alpha=0.995)
    x      = float.(x0)
    v      = f(x[1], x[2])
    T      = float(T0)
    x_best = copy(x)
    f_best = v
    v_prev = v

    for _ in 1:N
        # Korak 3: generiši susjedno rješenje
        dx    = delta_x .* (2 .* rand(2) .- 1)
        x_new = clamp.(x .+ dx,
                       [opseg[1][1], opseg[2][1]],
                       [opseg[1][2], opseg[2][2]])
        v_new = f(x_new[1], x_new[2])
        dv    = v_new - v

        # Metropolisov kriterij
        if dv < 0 || rand() < exp(-abs(dv) / T)
            x, v = copy(x_new), v_new
            if v < f_best
                x_best, f_best = copy(x), v
            end
        end

        # Korak 5: hlađenje
        T *= alpha

        # Korak 4: kriterij zaustavljanja
        if abs(v - v_prev) < eps && T < 1e-6
            break
        end
        v_prev = v
    end

    return x_best, f_best
end

# =========================
# Testne funkcije
# =========================
paraboloid(x,y) = x^2 + y^2
rastrigin(x,y)  = 20 + (x^2 - 10*cos(2*pi*x)) + (y^2 - 10*cos(2*pi*y))
dropwave(x,y)   = -(1 + cos(12*sqrt(x^2+y^2))) / (0.5*(x^2+y^2) + 2)
holder(x,y)     = -abs(sin(x)*cos(y)*exp(abs(1 - sqrt(x^2+y^2)/pi)))

# =========================
# Problem trgovačkog putnika (TSP)
# =========================
"""
    tour_length(cities, tour)

Ukupna dužina zatvorene ture `tour` kroz gradove čije su koordinate
redovi matrice `cities`.
"""
function tour_length(cities, tour)
    n = length(tour)
    s = 0.0
    for i in 1:n-1
        s += norm(cities[tour[i],:] - cities[tour[i+1],:])
    end
    s += norm(cities[tour[end],:] - cities[tour[1],:])
    return s
end

"""
    tsp_sa(cities; T0=1000.0, alpha=0.999, max_iter=200_000)

Simulirano kaljenje za TSP sa 2-opt operatorom susjedstva (obrtanje
dijela ture). Vraća par `(best_tour, best_len)`.
"""
function tsp_sa(cities; T0=1000.0, alpha=0.999, max_iter=200_000)
    n = size(cities,1)
    tour = shuffle(1:n)           # početna slučajna tura
    best = copy(tour)
    best_len = tour_length(cities, tour)
    T = T0

    for k in 1:max_iter
        i, j = sort(rand(1:n,2))
        new_tour = copy(tour)
        new_tour[i:j] = reverse(new_tour[i:j])

        old_len = tour_length(cities, tour)
        new_len = tour_length(cities, new_tour)
        Δ = new_len - old_len

        if Δ < 0 || rand() < exp(-Δ/T)
            tour = new_tour
            if new_len < best_len
                best = copy(new_tour)
                best_len = new_len
            end
        end

        T *= alpha
    end

    return best, best_len
end

end # module
