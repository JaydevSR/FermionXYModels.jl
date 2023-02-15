# FermionXYModel.jl

[![Build status (Github Actions)](https://github.com/JaydevSR/FermionXYModels.jl/workflows/CI/badge.svg)](https://github.com/JaydevSR/FermionXYModels.jl/actions)
[![codecov.io](http://codecov.io/github/JaydevSR/FermionXYModels.jl/coverage.svg?branch=main)](http://codecov.io/github/JaydevSR/FermionXYModels.jl?branch=main)

_**Quantum XY Model by projecting the spins to spinless fermions.**_

**Reference**: Stéphan, J., Misguich, G., &amp; Pasquier, V. (2010). Rényi entropy of a line in two-dimensional Ising models. Physical Review B, 82(12). doi:10.1103/physrevb.82.125455

## Installation

```julia
using Pkg
Pkg.add("FermionXYModels.jl")
```

## Models

- `FermionXYChain`: The quantum XY model as a chain of fermions.

    **Arguments:**
    - `L::Int`: The chain length.
    - `J::Real`: The coupling constant.
    - `h::Real`: The magnetic field.
    - `gamma::Real`: The anisotropy constant.
    - `start::Symbol`: The initial state. Can be one of `:rand`, `:vacuum`, `:filled`
    - `parity::Int=-1`: Can be -1 or 1.

- `FermionXXChain`: The quantum XY model for $\gamma=0$.

- `FermionIsingChain`: The quantum XY model for $\gamma=1$.

## Iterators

- `FermionBasis`: An Iterator over the binary basis of fermion chains. It accepts the following arguments:
   
   **Arguments:**
   - `n_sites::Int`: The number of sites in fermion chain.
   - `states::Tuple=(-1, 1)`: A keyword argument taking the integer representation of filled and unfilled sites.



## Correlations and Probabilities
- `correlation_matrix`: Calculates the correlation matrix of the chain given by $G_{ij} = \langle a_i^\dagger a_j\rangle$. Has two methods, one takes a `FermionXYChain` as argument. Other takes the agruments:

    - `L::Int`: The chain length.
    - `J::Real`: The coupling constant.
    - `h::Real`: The magnetic field.
    - `gamma::Real`: The anisotropy constant.
    - `parity::Int=-1`: Can be -1 or 1.
    - `float_type::Type=Float64`: The floating point type (default is `Float64`, can be `BigFloat` for greater precision).

- `probability_matrix`: Calculates the probability matrix of the chain. The probability of the particular configuration is then given by $\det(P)$ where $P$ is the said matrix. Has two methods, one takes a `FermionXYChain` as argument. Other takes the agruments:

    - `sites::Vector{Int}`: The sites of the chain having value -1 for no fermion and 1 for a fermion.
    - `L::Int`: The chain length.
    - `J::Real`: The coupling constant.
    - `h::Real`: The magnetic field.
    - `gamma::Real`: The anisotropy constant.
    - `parity::Int=-1`: Can be -1 or 1.
    - `float_type::Type=Float64`: The floating point type (default is `Float64`, can be `BigFloat` for greater precision).

## Monte-Carlo Simulation
- `metropolis_update!(model::FermionXYChain)`: Generates a new configuration for the chain by performing single site updates using acceptance rate $A(P'|P) = \cfrac{P'}{P}$, where $P'$ is the probability of new configuration and $P$ is that of old configuration..
- `equilibrate!(model::FermionXYChain, steps::Int)`: Equilibrates the chain by performing $N$ updates.

**Note: The above methods for monte-carlo sampling are useless from my experience this is probably because single site updates can not capture the transformation from quantum spins to spinless fermions. If anyone knows more about this please inform me by opening an issue**
