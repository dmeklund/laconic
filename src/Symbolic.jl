module Symbolic
    using ...Laconic: MatrixType
    import LinearAlgebra
    using LinearAlgebra: I

    abstract type AbstractStatement{T <: Tuple} end
    abstract type AbstractExpression{T <: Tuple} end
    struct Equals{T1,T2} <: AbstractStatement{Tuple{T1,T2}}
        lhs::T1
        rhs::T2
    end
    parenthesize(io::IO, expr::AbstractExpression) = print(io, "(", expr, ")")
    LinearAlgebra.dot(x::AbstractExpression, y::AbstractExpression) = x * y
    Base.transpose(expr::AbstractExpression) = expr
    Base.zero(expr::AbstractExpression) = Numeric(0)

    struct Variable <: AbstractExpression{Tuple{}}
        label::String
    end
    Base.show(io::IO, var::Variable) = print(io, var.label)
    parenthesize(io::IO, var::Variable) = print(io, var)
    # some Julia magic to make broadcasting work with our Variable type
    Base.length(var::AbstractExpression) = 1
    Base.iterate(var::AbstractExpression) = (var, nothing)
    Base.iterate(var::AbstractExpression, ::Any) = nothing
    Base.IteratorSize(::Type{<:AbstractExpression}) = Base.HasShape{0}()
    Base.size(x::AbstractExpression) = ()

    struct Numeric{T <: Number} <: AbstractExpression{Tuple{T}}
        value::T
    end
    Base.conj(expr::Numeric) = Numeric(conj(expr.value))
    parenthesize(io::IO, other) = print(io, other)
    parenthesize(io::IO, number::Numeric) = print(io, number)
    parenthesize(io::IO, number::Number) = print(io, number)
    parenthesize(io::IO, number::Numeric{Complex{T}}) where T = begin
        print(io, "(", number, ")")
    end
    parenthesize(io::IO, number::Numeric{Complex{Bool}}) = print(io, number)
    Base.show(io::IO, number::Numeric) = print(io, number.value)
    Base.:(==)(first::Numeric, second::Numeric) = (first.value == second.value)
    Base.:(==)(first::Numeric, second::Number) = (first.value == second)
    Base.:(==)(first::Number, second::Numeric) = (first == second.value)
    Base.one(elem::AbstractExpression) = Numeric{Integer}(1)

    struct NAryAddition{T <: Tuple} <: AbstractExpression{T}
        elements::T
    end
    Base.show(io::IO, expr::NAryAddition) = join(io, expr.elements, "+")
    NAryAddition(elem1, elem2) = begin
        #if combineterms(elem1) == combineterms(-elem2)
        #    Numeric(0)
        if elem1 == 0
            elem2
        elseif elem2 == 0
            elem1
        else
            NAryAddition((elem1, elem2))
        end
    end
    NAryAddition(elements...) = NAryAddition(elements)
    Base.conj(expr::NAryAddition) = begin
        NAryAddition([conj(element) for element in expr.elements])
    end

    struct Division{T1, T2} <: AbstractExpression{Tuple{T1,T2}}
        numerator::T1
        denominator::T2
    end
    Base.conj(expr::Division{T1,T2}) where {T1, T2} = Division(conj(expr.numerator), conj(expr.denominator))
    Base.show(io::IO, expr::Division{T1,T2}) where {T1,T2} = begin
        parenthesize(io, expr.numerator)
        print(io, "/")
        parenthesize(io, expr.denominator)
    end


    struct Exponential{T} <: AbstractExpression{Tuple{T}}
        argument::T
    end
    Base.conj(expr::Exponential) = Exponential(conj(expr.argument))
    Base.:(==)(lhs::Exponential, rhs::Exponential) = lhs.argument == rhs.argument
    parenthesize(io::IO, expr::Exponential) = print(io, expr)
    Base.show(io::IO, expr::Exponential{T}) where T = print(io, "exp(", expr.argument, ")")
    Exponential(expr::Numeric) = begin
        if expr == 0
            Numeric(1)
        else
            invoke(Exponential, Tuple{AbstractExpression}, expr)
        end
    end
    Base.exp(expr::AbstractExpression) = Exponential(expr)

    struct Max{T1, T2} <: AbstractExpression{Tuple{T1, T2}}
        arg1::T1
        arg2::T2
    end
    Base.max(arg1::AbstractExpression, arg2::AbstractExpression) = Max(arg1, arg2)
    Base.max(arg1::AbstractExpression, arg2::Number) = Max(arg1, Numeric(arg2))
    Base.max(arg1::Number, arg2::AbstractExpression) = Max(Numeric(arg1), arg2)
    Base.show(io::IO, expr::Max) = print(io, "max(", expr.arg1, ", ", expr.arg2, ")")

    struct Product{T <: Tuple} <: AbstractExpression{T}
        elements::T
    end
    Base.:(==)(lhs::Product, rhs::Product) = begin
        sort!([lhs.elements...], by=repr) == sort!([rhs.elements...], by=repr)
    end
    Base.conj(expr::Product) = Product([conj(element) for element in expr.elements]...)
    Base.copy(expr::Product{T}) where T = Product{T}(expr.elements)
    Base.show(io::IO, expr::Product{T}) where T = begin
        for (index, element) in enumerate(expr.elements)
            parenthesize(io, element)
            if index != length(expr.elements)
                print(io, "*")
            end
        end
    end
    Product(expr1::AbstractExpression{T1}, expr2::AbstractExpression{T2}) where {T1,T2} = begin
        if expr1 == 0 || expr2 == 0
            Numeric(0)
        elseif expr1 == 1
            expr2
        elseif expr2 == 1
            expr1
        elseif expr1 == expr2
            Power(expr1, Numeric(2))
        else
            Product((expr1, expr2))
        end
    end
    Product(arg1::AbstractExpression{T}, more...) where {T} = begin
        Product(tuple(arg1, more...))
    end
    Product(expr1::Product{T1}, expr2::Product{T2}) where {T1, T2} = begin
        Product(tuple(expr1.elements..., expr2.elements...))
    end
    Product(expr1::AbstractExpression{T1}, expr2::Product{T2}) where {T1,T2} = begin
        Product(expr1, expr2.elements...)
    end
    Product(expr1::Product{T1}, expr2::AbstractExpression{T2}) where {T1,T2} = begin
        Product(expr1.elements..., expr2)
    end
    Product(expr1::Exponential{T1}, expr2::Exponential{T2}) where {T1,T2} = begin
        Exponential(expr1.argument + expr2.argument)
    end
    Product(expr1::Product{Tuple{Exponential{T1}, T2}}, expr2::Product{Tuple{Exponential{T3}, T4}}) where {T1, T2, T3, T4} = begin
        (expr1.elements[1] * expr2.elements[1]) * expr1.elements[2] * expr2.elements[2]
    end

    struct Power{T1 <: AbstractExpression, T2 <: AbstractExpression} <: AbstractExpression{Tuple{T1, T2}}
        x::T1
        y::T2
    end
    Base.:(^)(x::T1, y::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = Power(x, y)
    Base.:(^)(x::T, y::Number) where {T <: AbstractExpression} = Power(x, Numeric(y))
    Base.conj(expr::Power{T1,Numeric{Integer}}) where T1 = Power(conj(expr.x), expr.y)
    Base.show(io::IO, expr::Power) = begin
        parenthesize(io, expr.x)
        print(io, "^")
        parenthesize(io, expr.y)
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
    Negation(expr::Numeric{T}) where T = Numeric(-expr.value)
    Negation(expr::Negation) = expr.element
    Negation(expr::Product) = Product(-expr.elements[1], expr.elements[2:end]...)
    Base.:(==)(lhs::Negation, rhs::Negation) = lhs.element == rhs.element
    Base.conj(expr::Negation) = Negation(conj(expr.element))
    Base.show(io::IO, expr::Negation{T}) where T = print(io, "-", expr.element)
    struct Inversion{T} <: AbstractArray{AbstractExpression{Tuple{T}},2}
        element::Array{T,2}
    end

    struct Abs{T} <: AbstractExpression{Tuple{T}}
        argument::T
    end
    Base.abs(argument::T) where T <: AbstractExpression = Abs(argument)

    Base.:*(obj::AbstractExpression) = obj
    Base.:+(obj::AbstractExpression) = obj
    Base.:*(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = Product(first, second)
    Base.:*(first::Number, second::T) where {T <: AbstractExpression} = Numeric(first) * second
    Base.:*(first::T, second::Number) where {T <: AbstractExpression} = first * Numeric(second)
    Base.:/(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = Division(first, second)
    Base.:/(first::T, second::Number) where {T <: AbstractExpression} = first / Numeric(second)
    Base.:+(first::NAryAddition, second::NAryAddition) = NAryAddition(first.elements..., second.elements...)
    Base.:+(first::AbstractExpression, second::NAryAddition) = NAryAddition(first, second.elements...)
    Base.:+(first::NAryAddition, second::AbstractExpression) = NAryAddition(first.elements..., second)
    Base.:+(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = NAryAddition(first, second)
    Base.:+(first::Number, second::AbstractExpression) = Numeric(first) + second
    Base.:+(first::AbstractExpression, second::Number) = first + Numeric(second)
    Base.:-(first::Number, second::AbstractExpression) = Numeric(first) - second
    Base.:-(first::T1, second::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression} = first + -second
    Base.:-(object::T) where {T <: AbstractExpression} = Negation(object)
    Base.:-(first::AbstractExpression, second::Number) = first - Numeric(second)

    LinearAlgebra.inv(expr::Array{T,2}) where {T <: AbstractExpression} = Inversion{T}(expr)
    Base.size(expr::Inversion{T}) where T = size(expr.element)

    struct Cosine{T} <: AbstractExpression{Tuple{T}}
        argument::T
    end
    Base.cos(expr::T) where {T <: AbstractExpression} = Cosine(expr)
    Base.conj(expr::Cosine) = Cosine(conj(expr.argument))
    parenthesize(io::IO, expr::Cosine{T}) where T = print(io, expr)
    Base.copy(expr::Cosine{T}) where T = Cosine{T}(expr.argument)
    Base.show(io::IO, expr::Cosine{T}) where T = print(io, "cos(", expr.argument, ")")

    struct Sine{T} <: AbstractExpression{Tuple{T}}
        argument::T
    end
    Base.sin(expr::T) where {T <: AbstractExpression} = Sine(expr)
    Base.conj(expr::Sine) = Sine(conj(expr.argument))
    parenthesize(io::IO, expr::Sine{T}) where T = print(io, expr)
    Base.show(io::IO, expr::Sine{T}) where T = print(io, "sin(", expr.argument, ")")

    struct SSymbol <: AbstractExpression{Tuple{String}}
        name::String
    end
    Base.conj(expr::SSymbol) = expr # FIXME
    parenthesize(io::IO, expr::SSymbol) = print(io, expr)
    Base.show(io::IO, symbol::SSymbol) = print(io, symbol.name)

    struct Vector{T <: AbstractArray} <: AbstractExpression{Tuple{T}}
        value::T
    end

    # simplification rules
    function NAryAddition(val1::Power{Sine{T1},Numeric{T2}}, val2::Power{Cosine{T3},Numeric{T4}}) where {T1,T2,T3,T4}
        if val1.y == 2 && val2.y == 2 && val1.x.argument == val2.x.argument
            Numeric{Integer}(1)
        else
            NAryAddition((val1, val2))
        end
    end
    function NAryAddition(val1::Power{Cosine{T1},Numeric{T2}}, val2::Power{Sine{T3},Numeric{T4}}) where {T1,T2,T3,T4}
        val2 + val1
    end

    # function combineterms(expr::Negation{Product{T}}) where T
    #     combineterms(Product(-expr.element.elements[1], expr.element.elements[2:end]...))
    # end

    function combineterms(expr::Product{T}) where T
        elements_copy = [combineterms(element) for element in expr.elements]
        sort!(elements_copy, by=repr)
        for (ind1, elem1) in enumerate(elements_copy)
            if isa(elem1, Negation)
                elements_copy[ind1] = elem1.element
                return Negation(combineterms(Product(elements_copy...)))
            end
            if isa(elem1, Product)
                return combineterms(Product(elements_copy[1:ind1-1]..., elem1.arguments..., elements_copy[ind1+1:end]))
            end
            for ind2 in ind1+1:length(elements_copy)
                elem2 = elements_copy[ind2]
                if elem1 == elem2
                    elements_copy[ind1] = Power(elem1, Numeric(2))
                    deleteat!(elements_copy, ind2)
                    return combineterms(Product(elements_copy...))
                end
                if isa(elem1, Negation) && isa(elem2, Negation)
                    elements_copy[ind1] = elem1.element
                    elements_copy[ind2] = elem2.element
                    return combineterms(Product(elements_copy...))
                end
                if isa(elem1, Exponential) && isa(elem2, Exponential)
                    elements_copy[ind1] = Exponential(elem1.argument + elem2.argument)
                    deleteat!(elements_copy, ind2)
                    return combineterms(Product(elements_copy...))
                end
            end
        end
        # sort!(elements_copy)
        Product(elements_copy...)
        # expr
    end
    function combineterms(expr::NAryAddition{T}) where T
        elements_copy = [expr.elements...] |> Array{AbstractExpression,1}
        for (ind1, elem1) in enumerate(elements_copy)
            elements_copy[ind1] = combineterms(elem1)
        end
        for (ind1, elem1) in enumerate(elements_copy)
            for (ind2, elem2) in enumerate(elements_copy[ind1+1:end])
                if (elem1 == -elem2)
                    elements_copy[ind1] = Numeric(0)
                    elements_copy[ind1+ind2] = Numeric(0)
                end
            end
        end
        NAryAddition(elements_copy...)
    end
    combineterms(expr::Exponential) = Exponential(combineterms(expr.argument))
    combineterms(expr) = expr

    function is_unitary(matrix::MatrixType{T}) where {T}
        product = matrix * conj(transpose(matrix))
        combineterms.(product) == I
    end

    convertToFunction(input, var1::Variable, var2::Variable) = begin
        (x) -> convertToFunction(convertToFunction(input, var1)(x[1]), var2)(x[2])
    end
    convertToFunction(input, var1::Variable, var2::Variable, var3::Variable) = begin
        (x) -> convertToFunction(convertToFunction(convertToFunction(input, var1)(x[1]), var2)(x[2]), var3)(x[3])
    end
    convertToFunction(input, variable::Variable) = x -> input
    convertToFunction(input::Variable, var::Variable) = begin
        if input == var
            x -> x
        else
            x -> input
        end
    end

    convertToFunction(expr::Power, var::Variable) = begin
        x -> convertToFunction(expr.x, var)(x) .^ convertToFunction(expr.y, var)(x)
    end
    convertToFunction(num::Numeric, var::Variable) = x -> num.value
    convertToFunction(expr::Product, var::Variable) = begin
        # need to use broadcast(*, ...) instead of prod(...) so that things like
        # "x * 3" work when x is a vector
        x -> broadcast(*, (convertToFunction(item, var)(x) for item in expr.elements)...)
    end
    convertToFunction(expr::NAryAddition, var::Variable) = begin
        # need to use broadcast(+, ...) instead of sum(...) so that things like
        # "x * 3" work when x is a vector
        x -> broadcast(+, (convertToFunction(item, var)(x) for item in expr.elements)...)
    end
    convertToFunction(expr::Sine, var::Variable) = begin
        x -> sin(convertToFunction(expr.argument, var)(x))
    end
    convertToFunction(expr::Division, var::Variable) = begin
        x -> convertToFunction(expr.numerator, var)(x) / convertToFunction(expr.denominator, var)(x)
    end
    convertToFunction(expr::Cosine, var::Variable) = begin
        x -> cos(convertToFunction(expr.argument, var)(x))
    end
    convertToFunction(expr::Negation, var::Variable) = begin
        x -> -convertToFunction(expr.element, var)(x)
    end
    convertToFunction(expr::Abs, var::Variable) = begin
        x -> abs(convertToFunction(expr.argument, var)(x))
    end
    convertToFunction(expr::Exponential, var::Variable) = begin
        x -> exp.(convertToFunction(expr.argument, var)(x))
    end
    convertToFunction(expr::Max, var::Variable) = begin
        x -> max(convertToFunction(expr.arg1, var)(x), convertToFunction(expr.arg2, var)(x))
    end

    evalexpr(expr, x::Variable, x0) = convertToFunction(expr, x)(x0)

    export AbstractStatement, AbstractExpression
    export Variable
    export Equals, NAryAddition
    export Numeric, SSymbol
    export Product, Exponential, Division, Power
    export Sine, Cosine
    export Negation, Abs
    export combineterms, is_unitary
    export convertToFunction, evalexpr
end
