# using Revise
# push!(LOAD_PATH, pwd())

"""
# module base.jl

- Julia version: 1.0
- Author: deklund
- Date: 2019-06-09

# Examples

```jldoctest
julia>
```
"""
module Laconic
    import LinearAlgebra: transpose, I
    import SparseArrays: sparse
    MatrixType{T} = AbstractArray{T, 2}
    VectorType{T} = AbstractArray{T, 1}

    struct Basis{T}
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
        fromBasisA = inv(toBasisA)
        basisA.basisToXform[basisB] = fromBasisA
        basisB
    end
    Basis{T}(name::String, ndims::Int64) where T = Basis{T}(name, ndims, nothing)
    Base.:(==)(basis1::Basis{T}, basis2::Basis{T}) where T = basis1.name == basis2.name && basis1.colnames == basis2.colnames && basis1.ndims == basis2.ndims

    function is_orthonormal(matrix::MatrixType{T}) where T
        product = transpose(matrix) * matrix
        product == I
    end

    struct Operator{T}
        name::String
        matrix::MatrixType{T}
        basis::Basis
        function Operator{T}(name::String, matrix::AbstractArray{T,2}, basis::Basis{T}) where T
            if size(matrix,1) != size(matrix,2)
                error("Matrix is not square: $(matrix)")
            elseif size(matrix,1) != basis.ndims
                error("Dimensions mismatch: $(size(matrix,1)) != $(basis.ndims)")
            else
                new{T}(name, matrix, basis)
            end
         end
    end
    transpose(op::Operator{T}, name::String) where T = Operator{T}(name, transpose(op.matrix), op.basis)
    transpose(op::Operator{T}) where T = transpose(op, "$(op.name).T")
    function Base.:+(op1::Operator{T}, op2::Operator{T}) where T
        if op1.basis != op2.basis
            error("Bases don't match: $(op1.basis), $(op2.basis)")
        end
        Operator{T}("$(op1.name)+$(op2.name)", op1.matrix+op2.matrix, op1.basis)
    end
    Operator{T}(name::String, op::Operator{T}) where T = Operator{T}(name, op.matrix, op.basis)
    function Base.:*(op::Operator{T}, scale::Number) where T
        Operator{T}("$(op.name)*$(scale)", op.matrix*scale, op.basis)
    end
    function Base.:/(op::Operator{T}, scale::Number) where {T}
        Operator{T}("$(op.name)/$(scale)", op.matrix/scale, op.basis)
    end
    function Base.:*(op1::Operator{T}, op2::Operator{T}) where T
        if op1.basis != op2.basis
            error("Bases don't match: $(op1.basis), $(op2.basis)")
        end
        Operator{T}("$(op1.name)*$(op2.name)", op1.matrix*op2.matrix, op1.basis)
    end
    function Base.:-(op1::Operator{T}, op2::Operator{T}) where T
        if op1.basis != op2.basis
            error("Bases don't match: $(op1.basis), $(op2.basis)")
        end
        Operator{T}("$(op1.name)-$(op2.name)", op1.matrix-op2.matrix, op1.basis)
    end
    function Base.:*(val::T2, op::Operator{T}) where {T, T2 <: Number}
        Operator{T}("$(val)*($op.name)", val*op.matrix, op.basis)
    end
    function Base.:(==)(op1::Operator{T}, op2::Operator{T}) where T
        op1.basis == op2.basis && op1.matrix ≈ op2.matrix
    end

    struct State{T}
        vector::VectorType{T}
        basis::Basis{T}
        function State{T}(vector::VectorType{T}, basis::Basis{T}) where T
            if basis.ndims == size(vector,1)
                new{T}(vector, basis)
            else
                error("Dimensions mismatch: $(basis.ndims) != $(size(vector,1))")
            end
        end
    end
    function convertToBasis(state::State{T}, basis::Basis{T}) where T
        if (state.basis == basis)
            state
        else
            xform = state.basis.basisToXform[basis]
            new_vec = xform * state.vector
            State{T}(new_vec, basis)
        end
    end

    function apply(operator::Operator{T}, state::State) where T
        if state.basis != operator.basis
            state = convertToBasis(state, operator.basis)
        end
        operator.matrix * state.vec
    end

    function commutator(op1::Operator{T}, op2::Operator{T}) where T
        op1*op2 - op2*op1
    end

    import LinearAlgebra: transpose
    import SparseArrays: sparse
    import Laconic: Basis, MatrixType, Operator

    struct Spin
        spin::Rational
        function Spin(spin::Rational)
            if spin.den in (1,2) && spin >= 0
                new(spin)
            else
                error("Invalid spin")
            end
        end
    end
    dicke(s::Spin) = Basis{ComplexF64}("Dicke", (2*s.spin).num+1)
    function splus(s::Spin)
        Is = 1:2*s.spin
        Js = 2:2*s.spin+1
        Ms = s.spin-1:-1:-s.spin
        Vs = sqrt.(s.spin*(s.spin+1) .- Ms.*(Ms.+1)) .+ 0.0im
        basis = dicke(s)
        Operator{ComplexF64}("splus", sparse(Is, Js, Vs, 2*s.spin+1, 2*s.spin+1), basis)
    end
    sminus(s::Spin) = transpose(splus(s), "sminus")
    sx(s::Spin) = Operator{ComplexF64}("sx", (splus(s) + sminus(s)) / 2.0)
    sy(s::Spin) = Operator{ComplexF64}("sy", (splus(s) - sminus(s)) / 2.0im)
    sz(s::Spin) = Operator{ComplexF64}(
        "sz",
        sparse(1:2*s.spin+1, 1:2*s.spin+1, s.spin:-1:-s.spin |> VectorType{ComplexF64}),
        dicke(s)
    )

end
