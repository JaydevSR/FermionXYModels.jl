module FermionXYModel

using LinearAlgebra
using ToeplitzMatrices

export FermionXYModel1D,
    FermionIsingModel1D,
    FermionXXModel1D,
    metropolis_update!,
    config_probability,
    correlation_matrix,
    probability_matrix

include("fermions.jl")
include("models.jl")
include("montecarlo.jl")

end # module
