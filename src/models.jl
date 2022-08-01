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
end

# -1 ≡ no fermion, +1 ≡ fermion
FermionXYChain(;L::Int, J::Real, h::Real, gamma::Real) = FermionXYChain{L}(
                                                                rand([-1, 1], L), L,
                                                                convert(Float64, J),
                                                                convert(Float64, h),
                                                                convert(Float64, gamma))

"""
The **XX model**, also known as the isotropic (γ=0) XY model, projected to spinless fermions with PBC.
"""
FermionXXChain(;L::Int, J::Real, h::Real) = FermionXYChain{L}(
                                rand([-1, 1], L), L,
                                convert(Float64, J),
                                convert(Float64, h),
                                0.0)
"""
The **Ising model**, i.e. the XY model with γ=1, projected to spinless fermions with PBC.
"""
FermionIsingChain(;L::Int, J::Real, h::Real) = FermionXYChain{L}(
                                    rand([-1, 1], L), L,
                                    convert(Float64, J),
                                    convert(Float64, h),
                                    1.0)

correlation_matrix(model::FermionXYChain; sign::Int=-1) = correlation_matrix(; 
                                                                        L=model.L, J=model.J,
                                                                        h=model.h, gamma=model.gamma,
                                                                        sign=sign)

function correlation_matrix(;L::Int, J::Float64, h::Float64, gamma::Float64, sign::Int=-1)
    return [G_nm(i, j; L=L, J=J, h=h, gamma=gamma, sign=sign) for i∈1:L, j∈1:L]
end

probability_matrix(model::FermionXYChain; sign=-1) = probability_matrix(model.sites;
                                                                    L=model.L, J=model.J,
                                                                    h=model.h, gamma=model.gamma,
                                                                    sign=sign)

@inbounds function probability_matrix(sites::Base.AbstractVecOrTuple; L::Int, J::Float64, h::Float64, gamma::Float64, sign::Int=-1)
    if length(sites) != L
        throw(ArgumentError("The number of sites is not equal to the model's length"))
    end
    return [(1/2)δ(i, j) - (1/2)*sites[i]*G_nm(i, j; L=L, J=J, h=h, gamma=gamma, sign=sign) for i∈1:L, j∈1:L]
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
