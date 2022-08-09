export FermionBasis

"""
    FermionBasis(n_sites; sites=(-1, 1))

Construct a binary fermion basis for a chain of fermions with `n_sites` in form of an interator. The default field is
`(-1, 1)` which can be optionally changed by specifying the keyword argument `states` as a `Tuple`. No allocations are made 
during construction. The basis can be materialized using [`collect`](@ref) (not recommended for long chains).
"""
struct FermionBasis
    n_sites::Int
    states::Tuple{Int,Int}
    function FermionBasis(n_sites::Int; states::Tuple=(-1, 1))
        if length(states) != 2
            throw(ArgumentError("length(states) != 2: Fermion basis requries two states. "))
        end
        if n_sites <= 0
            throw(ArgumentError("n_sites <= 0: Fermion basis requires n_sites > 0. "))
        end
        new(n_sites, states)
    end
end

function Base.iterate(b::FermionBasis, state::Int64=0)
    if state >= length(b)
        return nothing
    end
    indices = digits(state, base=2, pad=b.n_sites)
    return Tuple(b.states[i+1] for i in indices), state + 1
end

Base.length(b::FermionBasis) = 2^b.n_sites
