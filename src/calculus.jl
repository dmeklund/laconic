#=
calculus:
- Julia version: 1.2.0
- Author: david
- Date: 2019-12-26
=#
module Calculus
    using QuadGK
    using Laconic
    using Laconic.Symbolic

    struct DefiniteIntegral{T1, T2, T3} <: AbstractExpression{Tuple{T1,T2,T3}}
        variable::Variable
        startpoint::T1
        endpoint::T2
        integrand::T3
    end

    function evaluateIntegral(integral::DefiniteIntegral)
        integrandFunction = convertToFunction(integral.integrand, integral.variable)
        quadgk(integrandFunction, integral.startpoint, integral.endpoint)[1]
    end

    convertToFunction(input, variable::Variable) = x -> input
    convertToFunction(input::Variable, var::Variable) = begin
        if input == var
            x -> x
        else
            x -> input
        end
    end
    convertToFunction(integral::DefiniteIntegral, var::Variable) = begin
        if var == integral.variable
            error("Cannot make integral a function of the variable it's integrating over")
        end
        startpoint = convertToFunction(integral.startpoint, var)
        endpoint = convertToFunction(integral.endpoint, var)
        integrand = convertToFunction(integral.integrand, var)
        x -> evaluateIntegral(DefiniteIntegral(integral.variable, startpoint(x), endpoint(x), integrand(x)))
    end
    convertToFunction(expr::Power, var::Variable) where {T} = begin
        x -> convertToFunction(expr.x, var)(x) ^ convertToFunction(expr.y, var)(x)
    end
    convertToFunction(num::Numeric, var::Variable) = x -> num.value
    convertToFunction(expr::Product, var::Variable) = begin
        x -> prod(convertToFunction(item, var)(x) for item in expr.elements)
    end
    convertToFunction(expr::NAryAddition, var::Variable) = begin
        x -> sum(convertToFunction(item, var)(x) for item in expr.elements)
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

    evalexpr(expr, x::Variable, x0) = convertToFunction(expr, x)(x0)

    # TODO: move me somewhere else

    function positionfunc(basis::DiscreteMomentumBasis, n::Integer, x::Variable)
        sqrt(2/basis.a) * sin(n*pi*x/basis.a)
    end

    function positionfunc(basis::DiscretePositionBasis, n::Integer, x::Variable)
        result = 0
        mombasis = DiscreteMomentumBasis(basis.N, basis.a, basis.mass)
        delta = basis.a/(basis.N + 1)
        xj = delta * n
        for ind = 1:basis.N
            func = positionfunc(mombasis, ind, x)
            result = result + evalexpr(func, x, xj) * func
        end
        result
    end

    function integralidentity(r1::Variable, r2::Variable)
        t = Variable("t")
        2/sqrt(pi) * DefiniteIntegral(t, 0, Inf, Exponential(-t^2*(r1-r2)^2))
    end

    export Variable, DefiniteIntegral, convertToFunction, evaluateIntegral
    export positionfunc
    export evalexpr, integralidentity
end
