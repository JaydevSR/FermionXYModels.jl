"""
The **quantum XY model** projected to spinless fermions with **periodic boundary conditions** (PBC). The sites
consist of the elements from set {-1, 1}, where -1 represents no fermion and 1 represents a fermion at that site.
"""
mutable struct FermionXYModel1D{N}
    sites::Vector{Int}
    L::Int
    J::Float64
    h::Float64
    gamma::Float64
end

FermionXYModel1D(;L, J, h, gamma) = FermionXYModel1D{L}(fill(-1, L), L, J, h, gamma)

"""
The **XX model**, also known as the isotropic (γ=0) XY model, projected to spinless fermions with PBC.
"""
FermionXXModel1D(;L, J, h) = FermionXYModel1D{L}(fill(-1, L), L, J, h, 0)

"""
The **Ising model**, i.e. the XY model with γ=1, projected to spinless fermions with PBC.
"""
FermionIsingModel1D(;L, J, h) = FermionXYModel1D{L}(fill(-1, L), L, J, h, 1)

correlation_matrix(model::FermionXYModel1D; sign::Int=-1) = correlation_matrix(; 
                                                                        L=model.L, J=model.J,
                                                                        h=model.h, gamma=model.gamma,
                                                                        sign=sign)

probability_matrix(model::FermionXYModel1D; sign=-1) = probability_matrix(model.sites;
                                                                    L=model.L, J=model.J,
                                                                    h=model.h, gamma=model.gamma,
                                                                    sign=sign)

G_nm(model::FermionXYModel1D, n::Int, m::Int; sign::Int=-1) = G_nm(n, m;
                                                                L=model.L, J=model.J,
                                                                h=model.h, gamma=model.gamma,
                                                                sign=sign)
                                                  
function correlation_matrix(;L::Int, J::Float64, h::Float64, gamma::Float64, sign::Int=-1)
    return Toeplitz(
        [G_nm(n, 1; L=L, J=J, h=h, gamma=gamma, sign=sign) for n∈1:model.L],  # First column
        [G_nm(1, m; L=L, J=J, h=h, gamma=gamma, sign=sign) for m∈1:model.L])  # First row
end

function probability_matrix(sites::Base.AbstractVecOrTuple; L::Int, J::Float64, h::Float64, gamma::Float64, sign::Int=-1)
    if length(sites) != L
        throw(ArgumentError("The number of sites is not equal to the model's length"))
    end
    return Toeplitz(
        [(1/2)*δ(n,1) - (1/2)*model.sites[n]*G_nm(n, 1; L=L, J=J, h=h, gamma=gamma, sign=sign) for n∈1:model.L],
        [(1/2)*δ(1,m) - (1/2)*model.sites[1]*G_nm(n, 1; L=L, J=J, h=h, gamma=gamma, sign=sign) for m∈1:model.L])
end

function G_nm(n::Int, m::Int; L::Int, J::Float64, h::Float64, gamma::Float64, sign::Int=-1)
    g_n = 0
    for k in 1:L
        ϕ_k = (2k + (sign-1) // 2) // L  # redefinition of ϕ_k => ϕ_k / π
        ϵ_k = sqrt((J*cospi(ϕ_k) + h)^2 + (J*gamma*sinpi(ϕ_k))^2)
        cos_θ_k = (J * cospi(ϕ_k) + h) / ϵ_k
        sin_θ_k = J * gamma * sinpi(ϕ_k) / ϵ_k
        g_n += cos_θ_k*cospi((n-m)*ϕ_k) - sin_θ_k*sinpi((n-m)*ϕ_k)
    end
    return g_n/L
end

# Kronecker delta
δ(i, j) = isequal(i, j)
