#=
calculus:
- Julia version: 1.2.0
- Author: david
- Date: 2019-12-26
=#
module Calculus
    using Laconic
    using Laconic.Symbolic
    using Cubature
    using QuadGK

    struct DefiniteIntegral{T1, T2, T3} <: AbstractExpression{Tuple{T1,T2,T3}}
        variable::Variable
        startpoint::T1
        endpoint::T2
        integrand::T3
    end

    struct DefiniteIntegralN{N, T1, T2, T3} <: AbstractExpression{Tuple{T1,T2,T3}}
        vars::NTuple{N,Variable}
        startpoints::NTuple{N,T1}
        endpoints::NTuple{N,T2}
        integrand::T3
    end

    collapseintegrals(integral::DefiniteIntegral{T1, T2, DefiniteIntegral{T3,T4,T5}}) where {T1, T2, T3, T4, T5} = begin
        collapseintegrals(DefiniteIntegralN(
            (integral.variable, integral.integrand.variable),
            (integral.startpoint, integral.integrand.startpoint),
            (integral.endpoint, integral.integrand.endpoint),
            integral.integrand.integrand
        ))
    end
    collapseintegrals(integral::DefiniteIntegralN{N, T1, T2, DefiniteIntegral{T3,T4,T5}}) where {N, T1, T2, T3, T4, T5} = begin
        collapseintegrals(DefiniteIntegralN(
            (integral.vars..., integral.integrand.variable),
            (integral.startpoints..., integral.integrand.startpoint),
            (integral.endpoints..., integral.integrand.endpoint),
            integral.integrand.integrand
        ))
    end
    collapseintegrals(integral::DefiniteIntegral) = integral
    collapseintegrals(integral::DefiniteIntegralN) = integral

    Base.show(io::IO, integral::DefiniteIntegral) = print(io, "∫_$(integral.startpoint)^$(integral.endpoint) d$(integral.variable) $(integral.integrand)")

    evaluateintegral(integral::DefiniteIntegral) = begin
        println("Evaluating $(integral)")
        integrandFunction = convertToFunction(integral.integrand, integral.variable)
        quadgk(integrandFunction, integral.startpoint, integral.endpoint)[1]
    end
    evaluateintegral(integral::DefiniteIntegralN) = begin
        integrandfunction = convertToFunction(integral.integrand, integral.vars...)
        println("Evaluating $(integral.integrand)")
        pcubature(integrandfunction, integral.startpoints, integral.endpoints)[1]
    end


    Laconic.Symbolic.convertToFunction(integral::DefiniteIntegral, var::Variable) = begin
        if var == integral.variable
            error("Cannot make integral a function of the variable it's integrating over")
        end
        startpoint = convertToFunction(integral.startpoint, var)
        endpoint = convertToFunction(integral.endpoint, var)
        integrand = convertToFunction(integral.integrand, var)
        newintegral = x -> DefiniteIntegral(integral.variable, startpoint(x), endpoint(x), integrand(x))
        function (x)
            result = newintegral(x)
            # FIXME: sloppy way of handling this (evaluating the integral when
            # possible, otherwise just returning it as a definite integral)
            try
                evaluateintegral(result)
            catch err
                isa(err, MethodError) || rethrow(err)
                println("Integration failed, continuing")
                result
            end
        end
    end

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

    makefinite(expr::AbstractExpression) = expr

    function makefinite(integral::DefiniteIntegral)
        if integral.endpoint == Inf
            t = Variable("t0")
            f = convertToFunction(
                integral.integrand,
                integral.variable
            )
            if integral.startpoint == -Inf && integral.endpoint == Inf
                DefiniteIntegral(
                    t,
                    -1,
                    1,
                    f(t/(1-t)^2) * (1+t^2)/(1-t^2)^2
                )
            elseif isfinite(integral.startpoint) && integral.endpoint == Inf
                a = integral.startpoint
                DefiniteIntegral(
                    t,
                    0,
                    1,
                    f(a + t/(1-t)) / (1-t)^2
                )
            end
        elseif isfinite(integral.startpoint) && isfinite(integral.endpoint)
            integral
        else
            error("Case not currently handled: $(integral))")
        end
    end

    export Variable, DefiniteIntegral, DefiniteIntegralN, evaluateintegral
    export positionfunc
    export integralidentity, makefinite, collapseintegrals
    export convertToFunction
end
