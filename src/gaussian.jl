# based off pyquante2
module Gaussian
    function dist2(point1::NTuple{N,Float64}, point2::NTuple{N,Float64}) where N
        sum((point1 .- point2).^2)
    end

    struct PrimitiveGaussianBasisFunction{N}
        exponent::Float64
        powers::NTuple{N, Int64}
        origin::NTuple{N, Float64}
        normcoeff::Float64
    end

    function PrimitiveGaussianBasisFunction(
        exponent::Float64,
        powers::NTuple{N, Int64}=(0,0,0),
        origin::NTuple{N, Float64}=(0.,0.,0.),
    ) where N
        pgbf = PrimitiveGaussianBasisFunction(exponent, powers, origin, 1.0)
        norm = pgbfnorm(pgbf)
        PrimitiveGaussianBasisFunction(exponent, powers, origin, 1.0/norm)
    end

    pgbfnorm(pgbf) = sqrt(overlap(pgbf, pgbf))

    function overlap(
        a::PrimitiveGaussianBasisFunction,
        b::PrimitiveGaussianBasisFunction
    )
        gamma = a.exponent + b.exponent
        porigin = gaussianproductcenter(a.exponent, a.origin, b.exponent, b.origin)
        rab2 = dist2(a.origin, b.origin)
        pre = (pi/gamma)^1.5 * exp(-a.exponent*b.exponent*rab2/gamma)
        overlaps = (
            overlap1d(
                a.powers[ind],
                b.powers[ind],
                porigin[ind]-a.origin[ind],
                porigin[ind]-b.origin[ind],
                gamma)
            for ind=1:3
        )
        a.normcoeff * b.normcoeff * pre * prod(overlaps)
    end

    function gaussianproductcenter(expn1, center1, expn2, center2)
        (expn1 .* center1 .+ expn2 .* center2) ./ (expn1 + expn2)
    end

    function overlap1d(la::Int64, lb::Int64, ax::Float64, bx::Float64, gamma::Float64)
        total = 0
        for ind=0:div(la+lb, 2)
            total += binomial_prefactor(2ind, la, lb, ax, bx) * factorial2(2ind-1)/(2gamma)^ind
        end
        total
    end

    function amplitude(bf::PrimitiveGaussianBasisFunction, point)
        r2 = dist2(point, bf.origin)
        return bf.normcoeff * prod(r2 .^ bf.powers) * exp(-bf.exponent*r2)
    end

    function binomial_prefactor(s::Int64, ia::Int64, ib::Int64, xpa::Float64, xpb::Float64)
        total = 0
        for t=0:s
            if (s-ia) <= t <= ib
                total += binomial(ia, s-t)*binomial(ib, t)*xpa^(ia-s+t)*xpb^(ib-t)
            end
        end
        total
    end

    factorial2(n::Int64) = prod(n:-2:1)

    struct ContractedGaussianBasisFunction{N}
        powers::NTuple{3, Int64}
        origin::NTuple{3, Float64}
        normcoeff::Float64
        pgbfs::NTuple{N, PrimitiveGaussianBasisFunction}
        coeffs::NTuple{N, Float64}
    end

    function ContractedGaussianBasisFunction(
            powers,
            origin,
            pgbfs,
            coeffs
    )
        cgbf = ContractedGaussianBasisFunction(
            powers,
            origin,
            1.0,
            pgbfs,
            coeffs
        )
        norm = cgbfnorm(cgbf)
        ContractedGaussianBasisFunction(
            powers,
            origin,
            1.0/norm,
            pgbfs,
            coeffs
        )
    end

    cgbfnorm(cgbf::ContractedGaussianBasisFunction) = sqrt(overlap(cgbf, cgbf))

    overlap(a::ContractedGaussianBasisFunction, b::ContractedGaussianBasisFunction) = contract(overlap, a, b)

    primitives(a::ContractedGaussianBasisFunction) = zip(a.coeffs, a.pgbfs)

    function contract(func, a::ContractedGaussianBasisFunction, b::ContractedGaussianBasisFunction)
        result = 0
        for (coeffa, abf) in primitives(a)
            for (coeffb, bbf) in primitives(b)
                result += coeffa * coeffb * func(abf, bbf)
            end
        end
        return a.normcoeff * b.normcoeff * result
    end


    export PrimitiveGaussianBasisFunction, ContractedGaussianBasisFunction
    export amplitude
end
