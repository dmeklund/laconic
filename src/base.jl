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
module base
    import LinearAlgebra: transpose
    import SparseArrays: sparse
    MatrixType{T} = AbstractArray{T, 2}
    VectorType{T} = AbstractArray{T, 1}

    struct Basis
        name::String
        ndims::Int64
        basisToXform::AbstractDict{Basis,MatrixType{Any}}
        function Basis(name::String, ndims::Int64, basisToXform::AbstractDict{Basis,MatrixType{Any}}) where {T}
            new(name, ndims, basisToXform)
        end
    end
    Basis(name::String, ndims::Int64) = Basis(name, ndims, Dict{Basis,MatrixType{Any}}())
    Basis(name::String, basisA::Basis, toBasisA::MatrixType) = begin
        basisToXform = Dict{Basis,MatrixType}(basisA => toBasisA)
        basisB = Basis(name, basisA.ndims, basisToXform)
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
    function convertToBasis(state::State, basis::Basis)
        if (state.basis == basis)
            state
        else
            xform = state.basis.basisToXform[basis]
            new_vec = xform * state.vector |> VectorType
            State(new_vec, basis)
        end
    end

    function apply(operator::Operator, state::State)
        if state.basis != operator.basis
            state = convertToBasis(state, operator.basis)
        end
        operator.matrix * state.vec
    end

    function exercise2_1()
        increment(theta, psi) = [cos(theta/2), exp(im*psi)*sin(theta/2)] |> VectorType
        decrement(theta, psi) = [-exp(-im*psi)*sin(theta/2), cos(theta/2)] |> VectorType
        basisUpDown = Basis("UpDown", 2)
        for theta in (0, pi/4, pi)
            for psi in (0, pi/4, pi)
                inc = increment(theta, psi)
                dec = decrement(theta, psi)
                @assert dot(inc, dec) ≈ 0
                @assert dot(inc, inc) ≈ 1
                @assert dot(dec, dec) ≈ 1
                xform = hcat(inc, dec) |> MatrixType
                basisIncDec = Basis("IncDec", basisUpDown, xform)
                state_inc = State([1, 0] |> VectorType, basisIncDec)
                state_updown = convertToBasis(state_inc, basisUpDown)
                @assert state_updown.vector[1] ≈ inc[1] "$(state_updown.vector[1]) !≈ $(inc[1])"
                @assert state_updown.vector[2] ≈ inc[2] "$(state_updown.vector[2]) !≈ $(inc[2])"
            end
        end
    end

    module spin
        import LinearAlgebra: transpose
        import SparseArrays: sparse
        import base: Basis, MatrixType, Operator

        struct Spin
            spin::Rational
            Spin(spin::Rational) = spin.den in (1,2) && spin >= 0 ? new(spin) : error("Invalid spin")
        end
        dicke(s::Spin) = Basis("Dicke", 2*s.spin+1 |> Integer)
        function splus(s::Spin)
            Is = 1:2*s.spin
            Js = 2:2*s.spin+1
            Ms = s.spin-1:-1:-s.spin
            Vs = sqrt.(s.spin*(s.spin+1) .- Ms.*(Ms.+1))
            basis = dicke(s)
            Operator("splus", sparse(Is, Js, Vs, 2*s.spin+1, 2*s.spin+1) |> MatrixType, basis)
        end
        sminus(s::Spin) = transpose(splus(s), "sminus")
        sx(s::Spin) = Operator("sx", (splus(s) + sminus(s)) / 2)
        sy(s::Spin) = Operator("sy", (splus(s) - sminus(s)) / 2)
        sz(s::Spin) = Operator("sz", sparse(1:2*s.spin+1, 1:2*s.spin+1, s.spin:-1:-s.spin), dicke(s))
    end

end
