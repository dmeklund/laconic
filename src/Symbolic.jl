module Symbolic
    using Laconic: Basis, MatrixType, is_orthonormal
    import LinearAlgebra

    abstract type AbstractStatement{T <: Tuple} end
    abstract type AbstractExpression{T <: Tuple} end
    struct Equals{T1,T2} <: AbstractStatement{Tuple{T1,T2}}
        lhs::T1
        rhs::T2
    end
    Base.transpose(expr::AbstractExpression) = expr
    # struct BinaryAddition{T1, T2} <: AbstractExpression{Tuple{T1,T2}}
    #     first::T1
    #     second::T2
    # end
    struct NAryAddition{T <: Tuple} <: AbstractExpression{T}
        elements::T
    end
    NAryAddition(elements...) = NAryAddition(elements)
    struct Division{T1, T2} <: AbstractExpression{Tuple{T1,T2}}
        numerator::T1
        denominator::T2
    end
    struct Product{T <: Tuple} <: AbstractExpression{T}
        elements::T
    end
    Base.copy(expr::Product{T}) where T = Product{T}(expr.elements)
    struct Power{T1 <: AbstractExpression, T2 <: AbstractExpression} <: AbstractExpression{Tuple{T1, T2}}
        x::T1
        y::T2
    end
    Product(x1::T1, x2::T2) where {T1,T2} = begin
        if x1 == x2
            Power{T1,Numeric}(x1, Numeric{Integer}(2))
        else
            Product{Tuple{T1,T2}}((x1,x2))
        end
    end
    struct Negation{T} <: AbstractExpression{Tuple{T}}
        element::T
    end
    struct Inversion{T} <: AbstractArray{AbstractExpression{Tuple{T}},2}
        element::Array{T,2}
    end
    function Base.:*(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression}
        Product(first, second)
    end
    Base.:/(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = Division(first, second)
    Base.:+(first::NAryAddition, second::NAryAddition) = NAryAddition(first.elements..., second.elements...)
    function Base.:+(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression}
        NAryAddition(first, second)
    end
    Base.:-(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = T1 + -T2
    Base.:-(object::T) where {T <: AbstractExpression} = Negation(object)
    LinearAlgebra.inv(expr::Array{T,2}) where {T <: AbstractExpression} = Inversion{T}(expr)
    Base.size(expr::Inversion{T}) where T = size(expr.element)

    struct Cosine{T} <: AbstractExpression{Tuple{T}}
        argument::T
    end
    Base.copy(expr::Cosine{T}) where T = Cosine{T}(expr.argument)
    struct Sine{T} <: AbstractExpression{Tuple{T}}
        argument::T
    end
    struct Exponential{T} <: AbstractExpression{Tuple{T}}
        argument::T
    end
    struct Symbol <: AbstractExpression{Tuple{String}}
        name::String
    end
    struct Numeric{T <: Number} <: AbstractExpression{Tuple{T}}
        value::T
    end
    Base.:(==)(first::Numeric, second::Numeric) = (first.value == second.value)
    Base.:(==)(first::Numeric, second::Number) = (first.value == second)
    Base.:(==)(first::Number, second::Numeric) = (first == second.value)
    Base.one(elem::AbstractExpression) = Numeric{Integer}(1)
    struct Vector{T <: AbstractArray} <: AbstractExpression{Tuple{T}}
        value::T
    end
    function exercise2_1()
        state_up = Vector([1, 0])
        state_down = Vector([0, 1])
        theta = Symbol("ϑ")
        phi = Symbol("ϕ")
        increment(theta, phi) = [Cosine(theta/Numeric(2)), Exponential(Numeric(im)*phi)*Sine(theta/Numeric(2))]
        decrement(theta, phi) = [-Exponential(Numeric(-im)*phi)*Sine(theta/Numeric(2)), Cosine(theta/Numeric(2))]
        basisUpDown = Basis{AbstractExpression}("UpDown", 2, ("up", "down"))
        inc = increment(theta, phi)
        dec = decrement(theta, phi)
        xform = hcat(inc, dec)
        basisIncDec = Basis{AbstractExpression}(
            "IncDec",
            ("inc", "dec"),
            basisUpDown,
            xform |> MatrixType{AbstractExpression}
        )
        # @assert(is_orthonormal(xform))
        # print(xform)
        xform
        # Cosine(theta / Numeric(2)) * state_up + Exponential(Numeric(im)*phi) * Sine(theta / Numeric(2)) * state_down

    end

    # simplification rules
    function NAryAddition(val1::Power{Sine{T},Numeric}, val2::Power{Cosine{T},Numeric}) where T
        if val1.y == 2 && val2.y == 2 && val1.x.argument == val2.x.argument
            Numeric{Integer}(1)
        else
            NAryAddition((val1, val2))
        end
    end

    export AbstractStatement, AbstractExpression
    export Equals, BinaryAddition, NAryAddition
    export Numeric, Symbol
    export Product
    export Sine, Cosine
end
