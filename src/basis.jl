using Laconic.Symbolic

import LinearAlgebra

abstract type AbstractBasis end
abstract type AbstractOperator end
abstract type AbstractState end

hbar = 1

struct Basis{T} <: AbstractBasis
    name::String
    ndims::Int64
    colnames::Union{Nothing,Tuple{Vararg{AbstractString}}}
    basisToXform::AbstractDict{Basis{T},MatrixType{T}}
    function Basis{T}(
            name::String,
            ndims::Int64,
            colnames::Union{Nothing,Tuple{Vararg{AbstractString}}},
            basisToXform::AbstractDict{Basis{T},MatrixType{T}}
    ) where {T}
        if !isnothing(colnames) && length(colnames) != ndims
            error("Length mismatch: $(length(colname)) != $(ndims)")
        end
        new{T}(name, ndims, colnames, basisToXform)
    end
end

struct DiscreteBasis <: AbstractBasis
    indexToState::Function
    N::Integer
end

showInBasis(basis::Basis, showIn::Basis) = begin
    xform = basis.basisToXform[showIn]
    buffer = IOBuffer()
    # for eachrow in xform
    #
    # end
end

Basis{T}(name::String, ndims::Int64, colnames::Union{Nothing,Tuple{Vararg{AbstractString}}}) where {T} =
    Basis{T}(name, ndims, colnames, Dict{Basis{T},MatrixType{T}}())

function Basis{T}(
        name::String,
        colnames::Union{Nothing,Tuple{Vararg{AbstractString}}},
        basisA::Basis{T},
        toBasisA::MatrixType{T}
) where {T}
    basisToXform = Dict{Basis{T},MatrixType{T}}(basisA => toBasisA)
    basisB = Basis{T}(name, basisA.ndims, colnames, basisToXform)
    # fromBasisA = inv(toBasisA)
    # basisA.basisToXform[basisB] = fromBasisA
    basisB
end
Basis{T}(name::String, ndims::Int64) where T = Basis{T}(name, ndims, nothing)
Base.:(==)(basis1::Basis{T}, basis2::Basis{T}) where T = basis1.name == basis2.name && basis1.colnames == basis2.colnames && basis1.ndims == basis2.ndims
function LinearAlgebra.kron(basis1::Basis{T}, basis2::Basis{T}) where T
    Basis{T}("$(basis1.name)⊗$(basis2.name)", basis1.ndims*basis2.ndims, nothing)
end

struct MomentumBasis <: AbstractBasis
    a::AbstractFloat
end

function getBasisState(basis::MomentumBasis, index::Integer)
    MomentumEigenstate(basis.a, index, 1.0)
end

struct MomentumEigenstate <: AbstractState
    a::AbstractFloat
    n::Integer
    coeff::AbstractFloat
end

struct MomentumSquaredOperator <: AbstractOperator
end

struct PositionOperator <: AbstractOperator
end

struct ScaledOperator <: AbstractOperator
    coeff::AbstractFloat
    op::AbstractOperator
end

function apply(operator::MomentumSquaredOperator, state::MomentumEigenstate)
    MomentumEigenstate(
        state.a,
        state.n,
        state.coeff * state.n^2 * pi^2 * hbar^2 / state.a^2
    )
end

Base.:(*)(val::Number, op::AbstractOperator) = ScaledOperator(val, op)
Base.:(*)(op::AbstractOperator, val::Number) = val * op

function apply(operator::ScaledOperator, state::AbstractState)
    operator.coeff * apply(operator.op, state)
end

function createDiscreteBasis(basis::MomentumBasis, N::Integer)
    DiscreteBasis(n -> getBasisState(basis, n), N)
end

struct DiscretePositionBasis <: AbstractBasis
    N::Integer
    a::Real
    mass::Real
end

function symbolic(basis::DiscretePositionBasis, n::Int64, var::Variable)
    mombasis = DiscreteMomentumBasis(basis.N, basis.a, basis.mass)
    result = 0
    delta = basis.a/(basis.N+1)
    x_n = n * delta
    for ind=1:basis.N
        psiind = symbolic(mombasis, ind, var)
        func = convertToFunction(psiind, var)
        result = result + func(x_n) * psiind
    end
    sqrt(delta) * result
end

function xgrid(basis::DiscretePositionBasis)
    Array(1:basis.N) * basis.a / (basis.N + 1)
end

function createpos(basis::DiscretePositionBasis, x, sigma)
    vec = exp.(-((x .- xgrid(basis))/(2 * sigma)).^2)
    State(LinearAlgebra.normalize(vec), basis)
end

function getstate(basis::AbstractBasis, i::Integer)
    vec = zeros(ComplexF64, size(basis))
    vec[i] = 1.0
    State(vec, basis)
end

struct DiscreteMomentumBasis <: AbstractBasis
    N::Integer
    a::Real
    mass::Real
end

Base.size(basis::DiscreteMomentumBasis) = basis.N
Base.size(basis::DiscretePositionBasis) = basis.N

function symbolic(basis::DiscreteMomentumBasis, n::Int64, var::Variable)
    sqrt(2/basis.a) * Sine(n * pi * var / basis.a)
end

struct CombinedBasis{T <: Tuple} <: AbstractBasis
    bases::T
end

LinearAlgebra.kron(bases...) = CombinedBasis(bases)

function psix(basis::DiscreteMomentumBasis, n::Integer, x::Variable)
    sqrt(2/basis.a) * Sine(n*pi*x/basis.a)
end

export Basis, AbstractBasis, CombinedBasis
export AbstractOperator
export MomentumBasis, MomentumSquaredOperator, MomentumEigenstate
export PositionOperator
export DiscretePositionBasis, DiscreteMomentumBasis
export CombinedBasis
export createDiscreteBasis
export xgrid, psix, createpos
export getstate
export symbolic
