module Hydrogen
    using Laconic
    using Laconic.Calculus
    using Laconic.SystemM
    using Laconic.Symbolic

    struct HydrogenOrbital
        n::Int64
        l::Int64
        m::Int64
    end

    hbar = 1.0
    k_e = 1.0 # Coulomb's constant, 1/(4pi epsilon0)
    e = 1.0 # electron charge
    m_e = 1.0 # electron mass
    a0 = hbar^2 / (k_e * m_e * e^2) # Bohr radius

    function Laconic.symbolic(orb::HydrogenOrbital, r::Variable, theta::Variable, phi::Variable)
        normalization = sqrt((2/(orb.n*a0)^3) / factorial(orb.n + orb.l, orb.n - orb.l - 1) / (2*orb.n))
        ρ = 2*r / (orb.n * a0)
        psi = normalization * exp(-ρ/2) * ρ^orb.l * laguerre(2*orb.l+1, orb.n-orb.l-1, ρ) * harmonic(orb.l, orb.m, theta, phi)
        psi
    end

    function laguerre(alpha, n::Int64, x)
        if n < 0
            error("Invalid identifier: $(n)")
        elseif n == 0
            1
        elseif n == 1
            1 + alpha - x
        else
            k = n -1
            ((2*k + 1 + alpha - x) * laguerre(alpha, k, x) - (k + alpha) * laguerre(alpha, k-1, x)) / (k + 1)
        end
    end

    function harmonic(l::Int64, m::Int64, theta, phi)
        (-1)^m * sqrt((2*l+1)/4pi / factorial(l + m, l - m)) * legendre(l, m, cos(theta)) * exp(im*m*phi)
    end

    function legendre(l::Int64, m::Int64, x::Variable)
        1/(2^l * factorial(l)) * (1 - x^2)^(m/2) * deriv((x^2 - 1)^l, x, l + m)
    end

    function legendre(l::Int64, m::Int64, expr::AbstractExpression)
        x = Variable("x")
        result = legendre(l, m, x)
        convertToFunction(result, x)(expr)
    end

    export HydrogenOrbital, laguerre, harmonic, legendre
end