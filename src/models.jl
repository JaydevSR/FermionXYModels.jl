export FermionXYChain,
    FermionIsingChain,
    FermionXXChain,
    correlation_matrix,
    probability_matrix

"""
The **quantum XY model** projected to spinless fermions with **periodic boundary conditions** (PBC). The sites
consist of the elements from set {-1, 1}, where -1 represents no fermion and 1 represents a fermion at that site.
"""
mutable struct FermionXYChain{N}
    sites::Vector{Int}
    L::Int
    J::Float64
    h::Float64
    gamma::Float64
    corr::Matrix{Float64}
    sign::Int
    function FermionXYChain(; L::Int, J::Real, h::Real, gamma::Real, start::Symbol=:rand, sign::Int=-1)
        if L <= 0
            throw(ArgumentError("FermionXYChain requires L > 0."))
        end
        if !(-1 <= gamma <= 1)
            throw(ArgumentError("FermionXYChain requires -1 <= gamma <= 1."))
        end
        if J == 0
            throw(ArgumentError("FermionXYChain requires J to be non-zero."))
        end
        if sign != -1 && sign != 1
            throw(ArgumentError("FermionXYChain requires sign to be -1 or 1."))
        end
        J = convert(Float64, J)
        h = convert(Float64, h)
        gamma = convert(Float64, gamma)
        if start == :rand
            sites = rand([-1, 1], L)
        elseif start == :vacuum
            sites = ones(L)
        elseif start == :filled
            sites = fill(-1, L)
        else
            throw(ArgumentError("FermionXYChain allows start = :rand | :vacuum | :filled. "))
        end
        corr = correlation_matrix(; L=L, J=J, h=h, gamma=gamma, sign=sign)
        new{L}(sites, L, J, h, gamma, corr, sign)
    end
end

"""
The **XX model**, also known as the isotropic (γ=0) XY model, projected to spinless fermions with PBC.
"""
FermionXXChain(; L::Int, J::Real, h::Real, start::Symbol=:rand) = FermionXYChain(; L=L, J=J, h=h, gamma=0, start=start)

"""
The **Ising model**, i.e. the XY model with γ=1, projected to spinless fermions with PBC.
"""
FermionIsingChain(; L::Int, J::Real, h::Real, start::Symbol=:rand) = FermionXYChain(; L=L, J=J, h=h, gamma=1, start=start)

"""
    Returns the correlation matrix of the XY chain.
"""
correlation_matrix(model::FermionXYChain) = model.corr

"""
    Returns the probability matrix of the XY chain.
"""
probability_matrix(model::FermionXYChain) = probability_matrix(model.sites;
    L=model.L, J=model.J,
    h=model.h, gamma=model.gamma,
    sign=sign, corr=model.corr)

# correlation matrix of the fermion chain
function correlation_matrix(; L::Int, J::Float64, h::Float64, gamma::Float64, sign::Int=-1)
    return [G_nm(i, j; L=L, J=J, h=h, gamma=gamma, sign=sign) for i ∈ 1:L, j ∈ 1:L]
end

# probability matrix of the fermion chain
@inbounds function probability_matrix(sites::Base.AbstractVecOrTuple;
    L::Int,
    J::Float64,
    h::Float64,
    gamma::Float64,
    sign::Int=-1,
    corr::Matrix{Float64}=correlation_matrix(;
        L=L, J=J, h=h, gamma=gamma, sign=sign))

    if length(sites) != L
        throw(ArgumentError("The number of sites is not equal to the model's length"))
    end
    return [(1 / 2)δ(i, j) - (1 / 2) * sites[i] * corr[i, j] for i ∈ 1:L, j ∈ 1:L]
end

# elements of correlation matrix of the fermion chain
function G_nm(n::Int, m::Int; L::Int, J::Float64, h::Float64, gamma::Float64, sign::Int=-1)
    g_n = 0
    for k in 1:L
        ϕ_k = (2k + (sign - 1) // 2) // L  # redefinition of ϕ_k => ϕ_k / π
        ϵ_k = hypot(J * cospi(ϕ_k) + h, J * gamma * sinpi(ϕ_k))
        cos_θ_k = (J * cospi(ϕ_k) + h) / ϵ_k
        sin_θ_k = J * gamma * sinpi(ϕ_k) / ϵ_k
        g_n += cos_θ_k * cospi((n - m) * ϕ_k) - sin_θ_k * sinpi((n - m) * ϕ_k)
    end
    return g_n / L
end

# Kronecker delta
δ(i, j) = isequal(i, j)
