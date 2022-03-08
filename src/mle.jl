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

function PositivePreservingIterator(data::Matrix{T}, steps::S; dim::S) where {T<:Real, S<:Integer}
    π̂s = [Matrix{Complex{T}}(undef, dim, dim) for _ in 1:size(data, 2)]
    @sync for j in 1:size(data, 2)
        Threads.@spawn 𝛑̂!(π̂s[j], data[:, j]..., dim=dim)
    end

    ρ = glorot_uniform(Complex{T}, dim)
    ρ ./= tr(ρ)

    return PositivePreservingIterator{S, T}(dim, π̂s, ρ, steps)
end

"""
``\\sum_j \\frac{\\hat{\\pi_j}}{p_j}``
"""
function frac_π_p(ppit::PositivePreservingIterator{S, T}) where {S, T}
    sum_frac_πⱼ_pⱼ = zeros(Complex{T}, ppit.dim, ppit.dim)

    for π̂ in ppit.π̂s
        sum_frac_πⱼ_pⱼ .+= π̂ ./ tr_mul(π̂, ppit.ρ)
    end

    return sum_frac_πⱼ_pⱼ
end


function next!(ppit::PositivePreservingIterator)
    𝐫 = frac_π_p(ppit)

    ppit.ρ .= 𝐫 * ppit.ρ * 𝐫
    ppit.ρ ./= tr(ppit.ρ)

    return ppit.ρ
end

function run!(ppit::PositivePreservingIterator)
    p = Progress(ppit.steps)
    for _ in 1:ppit.steps
        next!(ppit)
        ProgressMeter.next!(p)
    end
end


# #########
# # utils #
# #########

glorot_uniform(T::Type{<:Real}, rng::AbstractRNG, dim) = (rand(rng, T, dim, dim) .- 0.5) .* sqrt(24 / 2dim)

glorot_uniform(::Type{Complex{U}}, rng, dim) where {U} = glorot_uniform(U, rng, dim) + im * glorot_uniform(U, rng, dim)

glorot_uniform(T, dim) = glorot_uniform(T, Random.GLOBAL_RNG, dim)

tr_mul(𝐚, 𝐛) = sum(𝐚[i, :]' * 𝐛[:, i] for i in 1:size(𝐚, 1))
