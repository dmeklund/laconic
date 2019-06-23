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
    import LinearAlgebra: transpose
    import SparseArrays: sparse
    MatrixType{T} = AbstractArray{T, 2}
    VectorType{T} = AbstractArray{T, 1}

    struct Basis{T}
        name::String
        ndims::Int64
        basisToXform::AbstractDict{Basis,MatrixType{T}}
        function Basis{T}(
                name::String,
                ndims::Int64,
                basisToXform::AbstractDict{Basis,MatrixType{T}}
        ) where {T}
            new{T}(name, ndims, basisToXform)
        end
    end

    Basis{T}(name::String, ndims::Int64) where {T} =
        Basis{T}(name, ndims, Dict{Basis,MatrixType{T}}())

    function Basis{T}(
            name::String,
            basisA::Basis{T},
            toBasisA::MatrixType{T}
    ) where {T}
        basisToXform = Dict{Basis,MatrixType{T}}(basisA => toBasisA)
        basisB = Basis{T}(name, basisA.ndims, basisToXform)
        fromBasisA = inv(toBasisA)
        basisA.basisToXform[basisB] = fromBasisA
        basisB
    end

    struct Operator
        name::String
        matrix::MatrixType
        basis::Basis
        function Operator(name::String, matrix::MatrixType, basis::Basis)
            if size(matrix,1) != size(matrix,2)
                error("Matrix is not square: $(matrix)")
            elseif size(matrix,1) != basis.ndims
                error("Dimensions mismatch: $(size(matrix,1)) != $(basis.ndims)")
            else
                new(name, matrix, basis)
            end
         end
    end
    transpose(op::Operator, name::String) = Operator(name, transpose(op.matrix), op.basis)
    transpose(op::Operator) = transpose(op, "$(op.name).T")
    function Base.:+(op1::Operator, op2::Operator)
        if op1.basis != op2.basis
            error("Bases don't match: $(op1.basis), $(op2.basis)")
        end
        Operator("$(op1.name)+$(op2.name)", op1.matrix+op2.matrix, op1.basis)
    end
    Operator(name::String, op::Operator) = Operator(name, op.matrix, op.basis)
    function Base.:*(op::Operator, scale::Number) Operator("$(op.name)*$(scale)", op.matrix*scale, op.basis)
    end
    function Base.:/(op::Operator, scale::Number)
        Operator("$(op.name)/$(scale)", op.matrix/scale, op.basis)
    end

    struct State
        vector::VectorType
        basis::Basis
        function State(vector::VectorType, basis::Basis)
            if basis.ndims == size(vector,1)
                new(vector, basis)
            else
                error("Dimensions mismatch: $(basis.ndims) != $(size(vector,1))")
            end
        end
    end
    function convertToBasis(state::State, basis::Basis{T}) where T
        if (state.basis == basis)
            state
        else
            xform = state.basis.basisToXform[basis]
            new_vec = xform * state.vector |> VectorType{T}
            State(new_vec, basis)
        end
    end

    function apply(operator::Operator, state::State)
        if state.basis != operator.basis
            state = convertToBasis(state, operator.basis)
        end
        operator.matrix * state.vec
    end

    module spin
        import LinearAlgebra: transpose
        import SparseArrays: sparse
        import Laconic: Basis, MatrixType, Operator

        struct Spin
            spin::Rational
            Spin(spin::Rational) = spin.den in (1,2) && spin >= 0 ? new(spin) : error("Invalid spin")
        end
        dicke(s::Spin) = Basis{Rational}("Dicke", 2*s.spin+1 |> Integer)
        function splus(s::Spin)
            Is = 1:2*s.spin
            Js = 2:2*s.spin+1
            Ms = s.spin-1:-1:-s.spin
            Vs = sqrt.(s.spin*(s.spin+1) .- Ms.*(Ms.+1))
            basis = dicke(s)
            Operator("splus", sparse(Is, Js, Vs, 2*s.spin+1, 2*s.spin+1) |> MatrixType{Rational}, basis)
        end
        sminus(s::Spin) = transpose(splus(s), "sminus")
        sx(s::Spin) = Operator("sx", (splus(s) + sminus(s)) / 2)
        sy(s::Spin) = Operator("sy", (splus(s) - sminus(s)) / 2)
        sz(s::Spin) = Operator("sz", sparse(1:2*s.spin+1, 1:2*s.spin+1, s.spin:-1:-s.spin), dicke(s))
    end

end
