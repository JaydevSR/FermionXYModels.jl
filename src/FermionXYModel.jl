module FermionXYModel

using LinearAlgebra

export FermionXYModel1D,
    FermionIsingModel1D,
    FermionXXModel1D,
    metropolis_update!,
    config_probability,
    correlation_matrix,
    probability_matrix

include("model.jl")
include("montecarlo.jl")

end # module
