# QuantumStateDistributions

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://foldfelis-QO.github.io/QuantumStateDistributions.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://foldfelis-QO.github.io/QuantumStateDistributions.jl/dev)
[![Build Status](https://github.com/foldfelis-QO/QuantumStateDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/foldfelis-QO/QuantumStateDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/foldfelis-QO/QuantumStateDistributions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/foldfelis-QO/QuantumStateDistributions.jl)

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia-repl
pkg> add QuantumStateDistributions
```

## Quick start

Construct a squeezed state via
[`QuantumStateBase.jl`](https://github.com/foldfelis-QO/QuantumStateBase.jl)
and declare its quantum state distribution.

```julia-repl
julia> using QuantumStateDistributions, QuantumStateBase

julia> r, θ, dim = 0.8, π/2, 100;

julia> ρ = SqueezedState(r, θ, Matrix, dim=dim)
100×100 Matrix{ComplexF64}:
       0.7477+0.0im          0.0+0.0im  …  -1.77607e-24-5.91946e-10im  0.0+0.0im
          0.0+0.0im          0.0+0.0im              0.0+0.0im          0.0+0.0im
 -2.14974e-17+0.351079im     0.0+0.0im      2.77945e-10-8.16923e-25im  0.0+0.0im
          0.0+0.0im          0.0+0.0im              0.0+0.0im          0.0+0.0im
             ⋮                          ⋱
          0.0+0.0im          0.0+0.0im              0.0+0.0im          0.0+0.0im
 -1.77607e-24+5.91946e-10im  0.0+0.0im      4.68637e-19+0.0im          0.0+0.0im
          0.0+0.0im          0.0+0.0im              0.0+0.0im          0.0+0.0im

julia> d = GaussianStateBHD(ρ);
```

Sample a point from the quantum state distribution
in intensity-to-measurement-phase quadrature coordinate
measured by balanced homodyne detection:

```julia-repl
julia> rand(d)
2-element Vector{Float64}:
 0.8420476666965236
 1.6008878775912423
```

Sample `n` points from the quantum state distribution:

```julia-repl
julia> n = 4096;

julia> rand(d, n)
2×4096 Matrix{Float64}:
  0.0018714   0.0034182   0.00403972  0.00780472  …   6.27393  6.27811    6.27884
 -0.706334   -1.16179    -0.195581    0.174201       -0.60763  0.853457  -0.217017
```

![](docs/src/assets/demo.png)
