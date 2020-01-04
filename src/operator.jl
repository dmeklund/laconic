struct Operator{T} <: AbstractOperator
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
function LinearAlgebra.kron(op1::Operator{T}, op2::Operator{T}) where T
    Operator{T}(
        "$(op1.name)⊗$(op2.name)",
        kron(op1.matrix, op2.matrix),
        kron(op1.basis, op2.basis)
    )
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
function LinearAlgebra.kron(state1::State{T}, state2::State{T}) where T
    State{T}(kron(state1.vector))
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

export apply
