"""
The **quantum XY model** projected to spinless fermions with **periodic boundary conditions** (PBC).
"""
mutable struct FermionXYModel1D{N}
    sites::Vector{Int}
    L::Int
    J::Float64
    h::Float64
    gamma::Float64
end

# -1 ≡ no fermion, +1 ≡ fermion
FermionXYModel1D(;L, J, h, gamma) = FermionXYModel1D{L}(fill(-1, L), L, J, h, gamma)

"""
The **XX model**, also known as the isotropic (γ=0) XY model, projected to spinless fermions with PBC.
"""
FermionXXModel1D(;L, J, h) = FermionXYModel1D{L}(rand([-1, 1], L), L, J, h, 0)

"""
The **Ising model**, i.e. the XY model with γ=1, projected to spinless fermions with PBC.
"""
FermionIsingModel1D(;L, J, h) = FermionXYModel1D{L}(rand([-1, 1], L), L, J, h, 1)

function correlation_matrix(model::FermionXYModel1D; sign=-1)
    return [G_nm(model, i, j; sign=sign) for i∈1:model.L, j∈1:model.L]
end

function correlation_matrix(L::Int, J::Float64, h::Float64, gamma::Float64, sign::Int=-1)
    return [G_nm(i, j; L=L, J=J, h=h, gamma=gamma, sign=sign) for i∈1:L, j∈1:L]
end

@inbounds function probability_matrix(model::FermionXYModel1D; sign=-1)
    return [(1/2)δ(i, j) - (1/2)*model.sites[i]*G_nm(model, i, j; sign=sign) for i∈1:model.L, j∈1:model.L]
end

function probability_matrix(sites::Base.AbstractVecOrTuple; L::Int, J::Float64, h::Float64, gamma::Float64, sign::Int=-1)
    if length(sites) != L
        throw(ArgumentError("The number of sites is not equal to the model's length"))
    end
    return [(1/2)δ(i, j) - (1/2)*sites[i]*G_nm(i, j; L=L, J=J, h=h, gamma=gamma, sign=sign) for i∈1:L, j∈1:L]
end

function G_nm(model::FermionXYModel1D, n::Int, m::Int; sign::Int=-1)
    g_n = 0
    for k in 1:model.L
        ϕ_k = (2k + (sign-1) // 2) // model.L  # redefinition of ϕ_k => ϕ_k / π
        ϵ_k = hypot(model.J*cospi(ϕ_k) + model.h, model.J*model.gamma*sinpi(ϕ_k))
        cos_θ_k = (model.J * cospi(ϕ_k) + model.h) / ϵ_k
        sin_θ_k = model.J * model.gamma * sinpi(ϕ_k) / ϵ_k
        g_n += cos_θ_k*cospi((n-m)*ϕ_k) - sin_θ_k*sinpi((n-m)*ϕ_k)
    end
    return g_n/model.L
end

function G_nm(n::Int, m::Int; L::Int, J::Float64, h::Float64, gamma::Float64, sign::Int=-1)
    g_n = 0
    for k in 1:L
        ϕ_k = (2k + (sign-1) // 2) // L  # redefinition of ϕ_k => ϕ_k / π
        ϵ_k = hypot(J*cospi(ϕ_k) + h, J*gamma*sinpi(ϕ_k))
        cos_θ_k = (J * cospi(ϕ_k) + h) / ϵ_k
        sin_θ_k = J * gamma * sinpi(ϕ_k) / ϵ_k
        g_n += cos_θ_k*cospi((n-m)*ϕ_k) - sin_θ_k*sinpi((n-m)*ϕ_k)
    end
    return g_n/L
end

# Kronecker delta
δ(i, j) = isequal(i, j)
