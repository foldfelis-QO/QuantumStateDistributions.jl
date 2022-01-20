struct GaussianStateBHD{T}<:ContinuousUnivariateDistribution
    state::T
end

#=
rand(::AbstractRNG, d::UnivariateDistribution)
sampler(d::Distribution)
logpdf(d::UnivariateDistribution, x::Real)
cdf(d::UnivariateDistribution, x::Real)
quantile(d::UnivariateDistribution, q::Real)
minimum(d::UnivariateDistribution)
maximum(d::UnivariateDistribution)
insupport(d::UnivariateDistribution, x::Real)

mean(d::UnivariateDistribution)
var(d::UnivariateDistribution)
modes(d::UnivariateDistribution)
mode(d::UnivariateDistribution)
skewness(d::UnivariateDistribution)
kurtosis(d::Distribution, ::Bool)
entropy(d::UnivariateDistribution, ::Real)
mgf(d::UnivariateDistribution, ::Any)
cf(d::UnivariateDistribution, ::Any)
=#

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

function mean_of_π̂ₓ²(θs::AbstractVector{<:Number}, ρ::AbstractMatrix)
    # ⟨πₓ²⟩ = ⟨ââ exp(-2im θ) + â†â† exp(2im θ) + ââ† + â†â⟩ / 4
    # ⟨πₓ²⟩ = (exp(-2im θ)⟨â²⟩ + exp(2im θ)⟨â†²⟩ + 1 + 2⟨ââ†⟩) / 4
    # here, ⟨ââ† + â†â⟩ = 1 + 2⟨â†â⟩ due to the commutation relation
    return (
        exp.(-2im*θs) .* mean_of_annihilate²(ρ) .+
        exp.(2im*θs) .* mean_of_create²(ρ) .+
        1 .+ 2mean_of_create_annihilate(ρ)
    ) ./ 4
end

function mean_of_π̂ₓ(θs::AbstractVector{<:Number}, ρ::AbstractMatrix)
    # ⟨πₓ⟩ = ⟨â exp(-im θ) + â† exp(im θ)⟩ / 2
    # ⟨πₓ⟩ = (exp(-im θ)⟨â⟩ + exp(im θ)⟨â†⟩) / 2
    return (
        exp.(-im*θs) .* mean_of_annihilate(ρ) .+
        exp.(im*θs) .* mean_of_create(ρ)
    ) ./ 2
end
