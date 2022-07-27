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
FermionXYModel1D(;L, J, h, gamma) = FermionXYModel1D{L}(rand([-1, 1], L), L, J, h, gamma)

"""
The **XX model**, also known as the isotropic (γ=0) XY model, projected to spinless fermions with PBC.
"""
FermionXXModel1D(;L, J, h) = FermionXYModel1D{L}(rand([-1, 1], L), L, J, h, 0)

"""
The **Ising model**, i.e. the XY model with γ=1, projected to spinless fermions with PBC.
"""
FermionIsingModel1D(;L, J, h) = FermionXYModel1D{L}(rand([-1, 1], L), L, J, h, 1)

function correlation_matrix(model::FermionXYModel1D)
    return [G_n(mode, i-j) for i∈1:model.L, j∈1:model.L]
end

@inbounds function probability_matrix(model::FermionXYModel1D)
    return [(1/2)δ(i, j) + (1/2)*model.sites[i]*G_n(model, i-j) for i∈1:model.L, j∈1:model.L]
end

#! h=1 gives NaN
#TODO: add different expression for h=1
function G_n(model::FermionXYModel1D, n::Int)
    g_n = 0
    for k in 1:model.L
        ϕ_k = 2 * k // model.L  # N = 1 for T = 0, redefinition of ϕ_k => ϕ_k / π
        ϵ_k = hypot(model.J*cospi(ϕ_k) - model.h, model.J*model.gamma*sinpi(ϕ_k))
        cos_θ_k = (model.J * cospi(ϕ_k) - model.h) / ϵ_k
        sin_θ_k = model.J * model.gamma * sinpi(ϕ_k) / ϵ_k
        g_n += cos_θ_k*cospi(ϕ_k*n) - sin_θ_k*sinpi(ϕ_k*n)
    end
    return g_n/model.L
end

# Kronecker delta
δ(i, j) = ==(i, j)