mutable struct FermionXYModel1D{N}
    sites::Vector{Int}
    L::Int
    J::Float64
    gamma::Float64
end

FermionXYModel1D(L, J, gamma) = FermionXYModel1D{L}(rand([-1, 1], L), L, J, gamma)

function correlation_matrix(model::FermionXYModel1D)
    G_corr = zeros(model.L, model.L)
end

function corr_g(model::FermionXYModel1D, h::Float64, n::Int)
    g_n = 0
    for k in eachindex(model.sites)
        ϕ_k = 2pi * k / model.L  # N = 1 for T = 0
        ϵ_k = hypot(model.J*cos(ϕ_k) - h, model.J*model.gamma*sin(ϕ_k))
        cos_θ_k = (model.J * cos(ϕ_k) - h) / ϵ_k
        sin_θ_k = model.J * model.gamma * sin(ϕ_k) / ϵ_k
        g_n += (cos_θ_k + im*sin_θ_k)*cis(ϕ_k*n)
    end
    return real(g_n)/model.L
end