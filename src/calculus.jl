#=
calculus:
- Julia version: 1.2.0
- Author: david
- Date: 2019-12-26
=#
module Calculus
    using QuadGK
    using Laconic.Symbolic

    struct DefiniteIntegral{T1, T2, T3}
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
    convertToFunction(expr::Power{Variable, T}, var::Variable) where {T} = begin
        x -> convertToFunction(expr.x, var)(x) ^ convertToFunction(expr.y, var)(x)
    end
    convertToFunction(num::Numeric, var::Variable) = x -> num.value
    convertToFunction(expr::Product, var::Variable) = begin
        x -> prod(convertToFunction(item, var)(x) for item in expr.elements)
    end

    export Variable, DefiniteIntegral, convertToFunction, evaluateIntegral
end