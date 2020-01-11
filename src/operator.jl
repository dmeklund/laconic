using SparseArrays
import FFTW

struct Operator{T1,T2} <: AbstractOperator where {T2 <: AbstractBasis}
    name::String
    matrix::MatrixType{T1}
    basis::T2
    # function Operator{T1,T2}(name::String, matrix::AbstractArray{T1,2}, basis::T2) where {T1,T2}
    #     if size(matrix,1) != size(matrix,2)
    #         error("Matrix is not square: $(matrix)")
    #     elseif size(matrix,1) != basis.ndims
    #         error("Dimensions mismatch: $(size(matrix,1)) != $(basis.ndims)")
    #     else
    #         new{T1,T2}(name, matrix, basis)
    #     end
    #  end
end
transpose(op::Operator{T}, name::String) where T = Operator{T}(name, transpose(op.matrix), op.basis)
transpose(op::Operator{T}) where T = transpose(op, "$(op.name).T")
function Base.:+(op1::Operator, op2::Operator)
    if op1.basis != op2.basis
        error("Bases don't match: $(op1.basis), $(op2.basis)")
    end
    Operator("$(op1.name)+$(op2.name)", op1.matrix+op2.matrix, op1.basis)
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
function Base.:*(val::T, op::Operator) where {T <: Number}
    Operator("$(val)*($op.name)", val*op.matrix, op.basis)
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
    vector::Vector{T}
    basis::DiscreteBasis
    # function State{T}(vector::Vector{T}, basis::DiscreteBasis) where T
    #     if basis.N == size(vector,1)
    #         new{T}(vector, basis)
    #     else
    #         error("Dimensions mismatch: $(basis.N) != $(size(vector,1))")
    #     end
    # end
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

function apply(op::AbstractOperator, state::State)
    result = zeros(ComplexF64, state.basis.N)
    for ind in 1:state.basis.N
        basis_state = state.basis.indexToState(ind)
        result += state.vector[ind] * apply(op, basis_state)
    end
    result
end

function commutator(op1::Operator{T}, op2::Operator{T}) where T
    op1*op2 - op2*op1
end

function positionoperator(basis::DiscretePositionBasis)
    matrix = spzeros(basis.N, basis.N)
    for row in 1:basis.N
        for col in 1:basis.N
            if row == col
                matrix[row,col] = basis.a / 2
            elseif (row - col) % 2 != 0
                matrix[row,col] = -8*basis.a*row*col/(pi^2 * (row^2 - col^2)^2)
            end
        end
    end
    Operator("position", matrix, basis)
end

function convertoperator(op::Operator{T,DiscretePositionBasis}, basis::DiscreteMomentumBasis) where T
    if op.basis.N != basis.N || !(op.basis.a ≈ basis.a) || !(op.basis.mass ≈ basis.mass)
        error("Bases not conformable")
    end
    matrix = FFTW.r2r(op.matrix, FFTW.RODFT00)/(2*(basis.N+1))
    Operator(op.name, matrix, basis)
end

function positionoperator(basis::DiscreteMomentumBasis)
    basispos = DiscretePositionBasis(basis.N, basis.a, basis.mass)
    pospos = positionoperator(basispos)
    convertoperator(pospos, basis)
end

function kineticenergyoperator(basis::DiscreteMomentumBasis)
    hbar = 1.0
    elements = Array(1:basis.N).^2 * pi^2 * hbar^2 / (2 * basis.mass * basis.a^2)
    matrix = spdiagm(0 => elements)
    Operator("kineticenergy", matrix, basis)
end

export apply
export State, Operator
export positionoperator, convertoperator
export kineticenergyoperator
