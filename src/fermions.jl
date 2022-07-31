export FermionBasis

"""
    FermionBasis(n_sites; sites=(-1, 1))

Construct a binary fermion basis for a chain of fermions with `n_sites` in form of an interator. The default field is
`(-1, 1)` which can be optionally changed by specifying the keyword argument `states` as a `Tuple`. No allocations are made 
during construction. The basis can be materialized using [`collect`](@ref) (not recommended for long chains).
"""
struct FermionBasis{T}
    n_sites::Int
    states::Tuple{T, T}
    function FermionBasis(n_sites::Int; states::Tuple{T, T}=(-1, 1)) where T
        if length(states) != 2
            throw(ArgumentError("Only 2 states are allowed for a binary fermion basis."))
        end
        new{T}(n_sites, states)
    end
end

function Base.iterate(b::FermionBasis, state::Int64=0)
    if state >= length(b)
        return nothing
    end
    indices = digits(state, base=2, pad=b.n_sites)
    return [b.states[i+1] for i in indices], state+1
end

Base.length(b::FermionBasis) = 2^b.n_sites
