# Simulirano kaljenje (Simulated Annealing)

Seminarski rad iz predmeta **Numerički algoritmi** — Elektrotehnički fakultet, Univerzitet u Sarajevu.

Implementacija metaheurističkog algoritma **simuliranog kaljenja** u programskom jeziku [Julia](https://julialang.org/), testirana na standardnim benchmark funkcijama i na problemu trgovačkog putnika (TSPLIB instanca *Berlin52*).

📄 **Cijeli rad:** [docs/Simulirano_Kaljenje.pdf](docs/Simulirano_Kaljenje.pdf)

---

## Sadržaj

- [O algoritmu](#o-algoritmu)
- [Struktura repozitorija](#struktura-repozitorija)
- [Pokretanje](#pokretanje)
- [Rezultati](#rezultati)
- [Literatura](#literatura)

## O algoritmu

Simulirano kaljenje (Kirkpatrick, Gelatt i Vecchi, 1983) zasniva se na analogiji s fizikalnim procesom kaljenja materijala: tvar se zagrije, a zatim sporo hladi kako bi dostigla stanje minimalne energije.

| Fizikalni sistem            | Optimizacijski problem      |
|-----------------------------|-----------------------------|
| Stanje sistema              | Kandidat rješenja           |
| Energija stanja             | Vrijednost funkcije cilja   |
| Temperatura                 | Upravljački parametar $T$   |
| Globalni energetski minimum | Globalni optimum            |

**Metropolisov kriterij.** Bolje rješenje se uvijek prihvata, a lošije s vjerovatnoćom koja opada s pogoršanjem i temperaturom:

$$
P(\text{prihvati } \mathbf{x}') =
\begin{cases}
1, & \Delta E < 0 \\
e^{-\Delta E / T}, & \text{inače}
\end{cases}
\qquad \Delta E = E(\mathbf{x}') - E(\mathbf{x})
$$

**Raspored hlađenja.** Koristi se geometrijski raspored $T_{k+1} = \alpha \cdot T_k$, $\alpha \in (0, 1)$.

**Koraci algoritma** (funkcija dvije promjenljive):

1. Energetska funkcija je sama funkcija cilja: $E(x_1, x_2) = f(x_1, x_2)$.
2. Odabere se početno rješenje $\mathbf{x}_0$ unutar domene.
3. Susjedno rješenje se generiše kao $\mathbf{x}' = \mathbf{x} + \Delta\mathbf{x}$, $\Delta\mathbf{x} \in (-\delta, \delta)^2$, i projektuje na opseg pretrage; prihvata se po Metropolisovom kriteriju.
4. Algoritam staje nakon $N$ iteracija ili kada promjena funkcije padne ispod $\varepsilon$ (uz dovoljno nisku temperaturu).
5. Temperatura se smanjuje i postupak se ponavlja od koraka 3.

## Struktura repozitorija

```
.
├── src/
│   └── SimulatedAnnealing.jl   # SA algoritam, testne funkcije, TSP (2-opt)
├── scripts/
│   ├── test_funkcije.jl        # SA vs. Nelder-Mead + contour plotovi
│   └── tsp_berlin52.jl         # SA na Berlin52 TSP instanci
├── data/
│   └── berlin52.txt            # koordinate 52 grada (TSPLIB)
├── docs/
│   ├── Simulirano_Kaljenje.pdf # seminarski rad
│   └── images/                 # slike iz rada
└── Project.toml                # Julia zavisnosti
```

### Glavna funkcija

```julia
SA(f, opseg, x0, delta_x, N, eps, T0; alpha=0.995) -> (x_best, f_best)
```

| Parametar | Opis |
|-----------|------|
| `f`       | funkcija cilja oblika `f(x1, x2)` |
| `opseg`   | opseg pretraživanja `[(xmin, xmax), (ymin, ymax)]` |
| `x0`      | početna tačka |
| `delta_x` | maksimalni pomak po svakoj koordinati |
| `N`       | maksimalan broj iteracija |
| `eps`     | minimalna promjena funkcije između dvije iteracije |
| `T0`      | početna temperatura |
| `alpha`   | faktor hlađenja (podrazumijevano `0.995`) |

Primjer:

```julia
include("src/SimulatedAnnealing.jl")
using .SimulatedAnnealing

x, fx = SA(holder, [(-10, 10), (-10, 10)], [9.0, -9.0], [0.5, 0.5], 50_000, 1e-8, 3.0)
```

## Pokretanje

Potrebna je Julia 1.9 ili novija. Iz korijena repozitorija:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

```bash
julia --project=. scripts/test_funkcije.jl
```

```bash
julia --project=. scripts/tsp_berlin52.jl
```

Skripte ispisuju rezultate u konzolu, a plotove i LaTeX tabelu spremaju u folder `results/`.

## Rezultati

### Benchmark funkcije

Algoritam je testiran na četiri dvodimenzionalne funkcije i upoređen s Nelder-Mead metodom (`Optim.jl`, `Fminbox(NelderMead())`) iz nasumične početne tačke.

| Funkcija | Formula | Globalni minimum |
|----------|---------|------------------|
| Paraboloid   | $x_1^2 + x_2^2$ | $f(0,0) = 0$ |
| Rastrigin    | $20 + x_1^2 - 10\cos(2\pi x_1) + x_2^2 - 10\cos(2\pi x_2)$ | $f(0,0) = 0$ |
| Drop-Wave    | $-\dfrac{1 + \cos\left(12\sqrt{x_1^2 + x_2^2}\right)}{0.5(x_1^2 + x_2^2) + 2}$ | $f(0,0) = -1$ |
| Holder Table | $-\left\lvert\sin x_1 \cos x_2 \exp\left(\left\lvert 1 - \frac{\sqrt{x_1^2+x_2^2}}{\pi}\right\rvert\right)\right\rvert$ | $f(\pm 8.055, \pm 9.6646) \approx -19.2085$ |

![Contour plotovi testnih funkcija](docs/images/testne_funkcije.png)

*Crveni kružić — globalni minimum, plavi križić — najbolje SA rješenje.*

Odstupanje od globalnog minimuma (najbolji rezultat po funkciji, iz rada):

| Funkcija     | SA                     | Nelder-Mead            |
|--------------|-----------------------:|-----------------------:|
| Paraboloid   | $1.3 \times 10^{-5}$   | $6.7 \times 10^{-10}$  |
| Rastrigin    | $16.9$                 | $43.2$                 |
| Drop-Wave    | $2.0 \times 10^{-3}$   | $0.814$                |
| Holder Table | $4.3 \times 10^{-5}$   | $11.41$                |

**Zaključak:** na glatkoj, konveksnoj funkciji Nelder-Mead je precizniji, ali na višejamnim funkcijama (Drop-Wave, Holder Table) SA pouzdano pronalazi globalni minimum dok Nelder-Mead zaglavi u lokalnom. Rastrigin ostaje težak i za SA uz ove parametre. Puna tabela za sve kombinacije $T_0 \in \{1, 3, 5\}$ i $\alpha \in \{0.995, 0.998\}$ nalazi se u radu.

### Problem trgovačkog putnika — Berlin52

Konfiguracija je permutacija gradova, a susjedno rješenje se dobija **2-opt** potezom (obrtanjem dijela ture). Parametri: $T_0 = 1000$, $\alpha = 0.999$, $200\,000$ iteracija, 30 pokretanja.

![Berlin52 tura](docs/images/tsp_berlin52.png)

| | Dužina ture |
|---|---:|
| Poznati optimum | 7542 |
| Najbolje SA rješenje (30 pokretanja) | 7733.23 |
| Relativna greška | ≈ 2.54 % |

U radu su opisane i primjene na problem ruksaka i optimizaciju neuronskih mreža.

## Literatura

1. S. Kirkpatrick, C. D. Gelatt, M. P. Vecchi, *Optimization by simulated annealing*, Science 220(4598), 671–680, 1983.
2. D. Delahaye, S. Chaimatanan, M. Mongeau, *Simulated annealing: From basics to applications*, Handbook of Metaheuristics, Springer, 2018.
3. W. H. Press et al., *Numerical Recipes: The Art of Scientific Computing*, 3rd ed., Cambridge University Press, 2007.
4. S. Lin, *Computer solutions of the traveling salesman problem*, Bell System Technical Journal 44(10), 2245–2269, 1965.
5. [TSPLIB — Berlin52](https://github.com/mastqe/tsplib/blob/master/berlin52.tsp)

---

**Autor:** Ahmed Spahić
