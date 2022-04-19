export GaussianStateBHD

struct QuantumStateBHD{T<:AbstractMatrix}<:ContinuousMultivariateDistribution
    ρ::T
end

"""
    GaussianStateBHD

A Gaussian state distribution in intensity-to-measurement-phase quadrature coordinate
measured by balanced homodyne detection.

## Example

```julia-repl
julia> using QuantumStateBase

julia> d = GaussianStateBHD(SqueezedState(0.8, π, Matrix, dim=100));

julia> points = rand(d, 4096);
```
"""
struct GaussianStateBHD{T<:AbstractMatrix}<:ContinuousMultivariateDistribution
    ρ::T
end

Base.length(d::GaussianStateBHD) = 2

sampler(d::GaussianStateBHD) = d

Base.eltype(::GaussianStateBHD{T}) where {T} = eltype(T)

function Distributions._rand!(rng::AbstractRNG, d::GaussianStateBHD, p::AbstractVector)
    T = real(eltype(d))

    p = real(p)
    p[1] = 2π * rand(rng, T)
    p[2] = mean(d, p[1]) + std(d, p[1]) * randn(rng, T)

    return p
end

function Distributions._rand!(rng::AbstractRNG, d::GaussianStateBHD, ps::AbstractMatrix)
    T = real(eltype(d))
    n = size(ps, 2)

    ps = real(ps)
    ps[1, :] .= sort!(2π .* rand(rng, T, n))
    ps[2, :] .= mean(d, view(ps, 1, :)) .+ std(d, view(ps, 1, :)) .* randn(rng, T, n)

    return ps
end

# Distributions._logpdf(d::GaussianStateBHD, x::AbstractArray)

mean(d::GaussianStateBHD, θs) = real(mean_of_π̂ₓ(d.ρ, θs))

var(d::GaussianStateBHD, θs) = real(mean_of_π̂ₓ²(d.ρ, θs) .- mean_of_π̂ₓ(d.ρ, θs).^2)

std(d::GaussianStateBHD, θs) = real(.√(mean_of_π̂ₓ²(d.ρ, θs) .- mean_of_π̂ₓ(d.ρ, θs).^2))

# entropy(d::GaussianStateBHD)

# cov(d::GaussianStateBHD)

# ##### for Gaussian state in intensity-to-measurement-phase quadrature coordinate #####

# π̂ₓ = (â exp(-im θ) + â† exp(im θ)) / 2

function mean_of_create(ρ::AbstractMatrix)
    # ⟨ψ|â†|ψ⟩ = tr(â† ρ)
    # Creation = [
    #     0   0   0   …
    #     √1  0   0   …
    #     0   √2  0   …
    #     0   0   √3  …
    #     ⋮   ⋮   ⋮
    # ]
    diag_1_of_ρ = diag(ρ, 1)
    dim = length(diag_1_of_ρ)
    T = eltype(ρ)

    return sum(.√(T.(1:dim)) .* diag_1_of_ρ)
end

function mean_of_create²(ρ::AbstractMatrix)
    # ⟨ψ|â†â†|ψ⟩ = tr((â†)² ρ)
    # Creation² = [
    #     0       0       0     …
    #     0       0       0     …
    #     √1√2    0       0     …
    #     0       √2√3    0     …
    #     0       0       √3√4  …
    #     ⋮       ⋮       ⋮
    # ]
    diag_2_of_ρ = diag(ρ, 2)
    dim = length(diag_2_of_ρ)
    T = eltype(ρ)

    return sum(.√(T.(1:dim)) .* .√(T.(2:(dim+1))) .* diag_2_of_ρ)
end

function mean_of_annihilate(ρ::AbstractMatrix)
    # ⟨ψ|â|ψ⟩ = tr(â ρ)
    # Annihilation = [
    #     0   √1  0   0   …
    #     0   0   √2  0   …
    #     0   0   0   √3  …
    #     ⋮   ⋮   ⋮   ⋮
    # ]
    diag_n1_of_ρ = diag(ρ, -1)
    dim = length(diag_n1_of_ρ)
    T = eltype(ρ)

    return sum(.√(T.(1:dim)) .* diag_n1_of_ρ)
end

function mean_of_annihilate²(ρ::AbstractMatrix)
    # ⟨ψ|ââ|ψ⟩ = tr((â)² ρ)
    # Annihilation² = [
    #    0        0       √1√2    0       0     …
    #    0        0       0       √2√3    0     …
    #    0        0       0       0       √3√4  …
    #    ⋮        ⋮        ⋮       ⋮       ⋮
    # ]
    diag_n2_of_ρ = diag(ρ, -2)
    dim = length(diag_n2_of_ρ)
    T = eltype(ρ)

    return sum(.√(T.(1:dim)) .* .√(T.(2:(dim+1))) .* diag_n2_of_ρ)
end

function mean_of_create_annihilate(ρ::AbstractMatrix)
    # ⟨ψ|â†â|ψ⟩ = tr(â†â ρ)
    # Creation * Annihilation = [
    #     0   0   0   0   …
    #     0   1   0   0   …
    #     0   0   2   0   …
    #     0   0   0   3   …
    #     ⋮   ⋮   ⋮   ⋮
    # ]
    diag_0_of_ρ = diag(ρ)
    dim = length(diag_0_of_ρ)
    T = eltype(ρ)

    return sum(T.(0:(dim-1)) .* diag_0_of_ρ)
end

function mean_of_π̂ₓ²(ρ::AbstractMatrix, θs)
    # ⟨πₓ²⟩ = ⟨ââ exp(-2im θ) + â†â† exp(2im θ) + ââ† + â†â⟩ / 4
    # ⟨πₓ²⟩ = (exp(-2im θ)⟨â²⟩ + exp(2im θ)⟨â†²⟩ + 1 + 2⟨ââ†⟩) / 4
    # here, ⟨ââ† + â†â⟩ = 1 + 2⟨â†â⟩ due to the commutation relation
    return (
        exp.(-2im*θs) .* mean_of_annihilate²(ρ) .+
        exp.(2im*θs) .* mean_of_create²(ρ) .+
        1 .+ 2mean_of_create_annihilate(ρ)
    ) ./ 4
end

function mean_of_π̂ₓ(ρ::AbstractMatrix, θs)
    # ⟨πₓ⟩ = ⟨â exp(-im θ) + â† exp(im θ)⟩ / 2
    # ⟨πₓ⟩ = (exp(-im θ)⟨â⟩ + exp(im θ)⟨â†⟩) / 2
    return (
        exp.(-im*θs) .* mean_of_annihilate(ρ) .+
        exp.(im*θs) .* mean_of_create(ρ)
    ) ./ 2
end
