module FermionXYModel

using LinearAlgebra

export FermionXYModel1D, metropolis_update!

include("model.jl")
include("montecarlo.jl")

end # module
