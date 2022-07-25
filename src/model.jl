mutable struct FermionXYModel1D{N}
    sites::Vector{Int}
    L::Int
    h::Float64
    J::Float64
    gamma::Float64
end

QuantumXYModel1D(L, h, J, gamma) = FermionXYModel1D{L}(rand([-1, 1], L), L, h, J, gamma)

function correlation_matrix(model::FermionXYModel1D)
    
end

function corr_g(model::FermionXYModel1D, n::Int)
    g_n = 0
    for k in eachindex(model.sites)
        ϕ_k = 2pi * k / model.L  # N = 1 for T = 0
        ϵ_k = sqrt((model.J * cos(ϕ_k) - h)^2 + model.J^2 * model.gamma^2 * sin(ϕ_k)^2)
        θ_k = acos((model.J * cos(ϕ_k) - h) / ϵ_k)
        g_n += cis(2n * pi * k / model.L + θ_k)
    end
    return g_n/model.L
end