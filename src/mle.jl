export
    PositivePreservingIterator,
    next!,
    run!

mutable struct PositivePreservingIterator{S, T}
    dim::S
    π̂s::Vector{Matrix{Complex{T}}}
    ρ::Matrix{Complex{T}}
    steps::S
end

function gen_π̂s(θs::AbstractVector{T}, xs::AbstractVector{T}, dim) where {T<:Real}
    n = length(θs)
    π̂s = [Matrix{Complex{T}}(undef, dim, dim) for _ in 1:n]

    return gen_π̂s!(π̂s, θs, xs, dim)
end

function gen_π̂s!(π̂s, θs::AbstractVector, xs::AbstractVector, dim)
    n = length(θs)

    @info "preprocessing..."
    p = Progress(n)
    Threads.@threads for i in 1:n
        𝛑̂!(π̂s[i], θs[i], xs[i], dim=dim)
        ProgressMeter.next!(p)
    end

    return π̂s
end

function PositivePreservingIterator(data::Matrix{T}, steps::S; dim::S) where {T<:Real, S<:Integer}
    π̂s = gen_π̂s(data[1, :], data[2, :], dim)
    ρ = glorot_uniform(Complex{T}, dim)
    ρ ./= tr(ρ)

    return PositivePreservingIterator{S, T}(dim, π̂s, ρ, steps)
end

"""
    frac_π_p(ppit::PositivePreservingIterator)

A positive preserving iterator.

``\\mathcal{R} = \\sum_i \\frac{\\hat{\\pi_i}}{p_i}``
"""
function frac_π_p(ppit::PositivePreservingIterator{S, T}) where {S, T}
    sum_frac_πᵢ_pᵢ = zeros(Complex{T}, ppit.dim, ppit.dim)

    for π̂ in ppit.π̂s
        sum_frac_πᵢ_pᵢ .+= π̂ ./ tr_mul(π̂, ppit.ρ)
    end

    return sum_frac_πᵢ_pᵢ
end

"""
    next!(ppit::PositivePreservingIterator)

Iterate one step with the magic positive preserving iterator.

``\\rho^{t+1} = \\mathcal{R} \\rho^t \\mathcal{R}``
"""
function next!(ppit::PositivePreservingIterator)
    𝐫 = frac_π_p(ppit)

    ppit.ρ .= 𝐫 * ppit.ρ * 𝐫
    ppit.ρ ./= tr(ppit.ρ)

    return ppit
end

"""
    run!(ppit::PositivePreservingIterator)

Iterate `n` step with the magic positive preserving iterator.
"""
function run!(ppit::PositivePreservingIterator)
    @info "estimating..."
    p = Progress(ppit.steps)
    for _ in 1:ppit.steps
        next!(ppit)
        ProgressMeter.next!(p)
    end

    return ppit
end

# #########
# # utils #
# #########

glorot_uniform(T::Type{<:Real}, rng::AbstractRNG, dim) = (rand(rng, T, dim, dim) .- 0.5) .* sqrt(24 / 2dim)

glorot_uniform(::Type{Complex{U}}, rng, dim) where {U} = glorot_uniform(U, rng, dim) + im * glorot_uniform(U, rng, dim)

glorot_uniform(T, dim) = glorot_uniform(T, Random.GLOBAL_RNG, dim)

tr_mul(𝐚, 𝐛) = sum(𝐚[i, :]' * 𝐛[:, i] for i in 1:size(𝐚, 1))
