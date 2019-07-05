module Symbolic
    export AbstractStatement, AbstractExpression, Equals, BinaryAddition, NAryAddition
    using Laconic: Basis, MatrixType
    import LinearAlgebra

    abstract type AbstractStatement{T <: Tuple} end
    abstract type AbstractExpression{T <: Tuple} end
    struct Equals{T1,T2} <: AbstractStatement{Tuple{T1,T2}}
        lhs::T1
        rhs::T2
    end
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
    struct Negation{T} <: AbstractExpression{Tuple{T}}
        element::T
    end
    struct Inversion{T} <: AbstractExpression{Tuple{T}}
        element::Array{T,2}
    end
    Base.:*(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = Product((first, second))
    Base.:/(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = Division(first, second)
    Base.:+(first::NAryAddition, second::NAryAddition) = NAryAddition(first.elements..., second.elements...)
    Base.:+(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = NAryAddition((first, second))
    Base.:-(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = T1 + -T2
    Base.:-(object::T) where {T <: AbstractExpression} = Negation(object)
    LinearAlgebra.inv(expr::Array{T,2}) where {T <: AbstractExpression} = Inversion{T}(expr)

    struct Cosine{T} <: AbstractExpression{Tuple{T}}
        argument::T
    end
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

        # Cosine(theta / Numeric(2)) * state_up + Exponential(Numeric(im)*phi) * Sine(theta / Numeric(2)) * state_down

    end
end
