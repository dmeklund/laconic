# based off pyquante2
module Gaussian
    using Laconic
    using Laconic.Symbolic
    using Distributed
    using SharedArrays
    using SparseArrays
    using SpecialFunctions

    function dist2(point1::NTuple{N,Float64}, point2::NTuple{N,Float64}) where N
        sum((point1 .- point2).^2)
    end

    struct PrimitiveGaussianBasisFunction{N}
        exponent::Float64
        powers::NTuple{N, Int64}
        origin::NTuple{N, Float64}
        normcoeff::Float64
    end

    function Laconic.symbolic(pgbf::PrimitiveGaussianBasisFunction{N}, x::NTuple{N,Variable}) where N
        pgbf.normcoeff * prod(x .^ pgbf.powers) * exp(-pgbf.exponent * sum((x .- pgbf.origin).^2))
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

    PrimitiveGaussianBasisFunction(
        exponent::Float64,
        N::Integer
    ) = PrimitiveGaussianBasisFunction(
        exponent, 
        ntuple(i->0, N),
        ntuple(i->0.0, N)
    )

    pgbfnorm(pgbf) = sqrt(overlap(pgbf, pgbf))

    function overlap(
        a::PrimitiveGaussianBasisFunction,
        b::PrimitiveGaussianBasisFunction
    )
        a.normcoeff * b.normcoeff * overlap(
            a.exponent,
            a.powers,
            a.origin,
            b.exponent,
            b.powers,
            b.origin
        )
    end

    function overlap(
            exponent1::Float64,
            powers1::NTuple{N,Int64},
            origin1::NTuple{N,Float64},
            exponent2::Float64,
            powers2::NTuple{N,Int64},
            origin2::NTuple{N,Float64},
    ) where N
        gamma = exponent1 + exponent2
        porigin = gaussianproductcenter(exponent1, origin1, exponent2, origin2)
        rab2 = dist2(origin1, origin2)
        pre = (pi/gamma)^1.5 * exp(-exponent1*exponent2*rab2/gamma)
        overlaps = (
            overlap1d(
                powers1[ind],
                powers2[ind],
                porigin[ind] - origin1[ind],
                porigin[ind] - origin2[ind],
                gamma)
            for ind=1:N
        )
        pre * prod(overlaps)
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
        delta = point .- bf.origin
        r2 = dist2(point, bf.origin)
        return bf.normcoeff * prod(delta .^ bf.powers) * exp(-bf.exponent*r2)
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

    struct ContractedGaussianBasisFunction{N, M}
        powers::NTuple{N, Int64}
        origin::NTuple{N, Float64}
        normcoeff::Float64
        pgbfs::NTuple{M, PrimitiveGaussianBasisFunction{N}}
        coeffs::NTuple{M, Float64}
    end

    function Laconic.symbolic(cgbf::ContractedGaussianBasisFunction{N}, vars::NTuple{N,Variable}) where N
        cgbf.normcoeff * sum(coeff * symbolic(pgbf, vars) for (coeff, pgbf) in zip(cgbf.coeffs, cgbf.pgbfs))
    end

    function ContractedGaussianBasisFunction(
            powers::NTuple{N,Int64},
            origin::NTuple{N,Float64},
            pgbfs::NTuple{M, PrimitiveGaussianBasisFunction},
            coeffs::NTuple{M, Float64}
    ) where {N, M}
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

    function ContractedGaussianBasisFunction(
            powers::NTuple{N,Int64},
            origin::NTuple{N,Float64},
            exponentsAndCoeffs::Tuple{Float64,Float64}...
    ) where {N}
        pgbfs = tuple(
            (PrimitiveGaussianBasisFunction(exponent, powers, origin)
             for (exponent, coeff) in exponentsAndCoeffs)...
         )
        coeffs = tuple((coeff for (exponent, coeff) in exponentsAndCoeffs)...)
        ContractedGaussianBasisFunction(powers, origin, pgbfs, coeffs)
    end

    cgbfnorm(cgbf::ContractedGaussianBasisFunction) = sqrt(overlap(cgbf, cgbf))

    overlap(a::ContractedGaussianBasisFunction, b::ContractedGaussianBasisFunction) = contract(overlap, a, b)

    primitives(a::ContractedGaussianBasisFunction{N}) where N = zip(a.coeffs, a.pgbfs)

    function contract(
        func,
        a::ContractedGaussianBasisFunction{N, M},
        b::ContractedGaussianBasisFunction{N, M}
    ) where {N, M}
        result = 0.
        for ind1=1:N
            for ind2=1:N
                result += a.coeffs[ind1] * b.coeffs[ind2] * func(a.pgbfs[ind1], b.pgbfs[ind2])::Float64
            end
        end
        # for (coeffa, abf) in primitives(a)
        #     for (coeffb, bbf) in primitives(b)
        #         result += coeffa * coeffb * func(abf, bbf)
        #     end
        # end
        return a.normcoeff * b.normcoeff * result
    end

    function contract(
            f,
            a::ContractedGaussianBasisFunction{N, M},
            b::ContractedGaussianBasisFunction{N, M},
            c::ContractedGaussianBasisFunction{N, M},
            d::ContractedGaussianBasisFunction{N, M}
    ) where {N,M}
        s = 0.
        for (ca, abf) in primitives(a)
            for (cb, bbf) in primitives(b)
                for (cc, cbf) in primitives(c)
                    for (cd, dbf) in primitives(d)
                        s += ca*cb*cc*cd*f(abf,bbf,cbf,dbf)
                    end
                end
            end
        end
        return a.normcoeff * b.normcoeff * c.normcoeff * d.normcoeff * s
    end


    function amplitude(
            bf::ContractedGaussianBasisFunction{N,M},
            point::NTuple{N, Float64}
    ) where {N, M}
        s = sum(c*amplitude(pbf, point) for (c, pbf) in primitives(bf))
        bf.normcoeff * s
    end

    function kinetic(
            pgbf1::PrimitiveGaussianBasisFunction{N},
            pgbf2::PrimitiveGaussianBasisFunction{N}
    ) where N
        overlap0 = overlap(
            pgbf1.exponent,
            pgbf1.powers,
            pgbf1.origin,
            pgbf2.exponent,
            pgbf2.powers,
            pgbf2.origin
        )
        term0 = pgbf2.exponent * (2 * sum(pgbf2.powers) + 3) * overlap0
        term1 = 0.
        term2 = 0.
        for ind=1:N
            term1 += -2 * (pgbf2.exponent^2) * overlap(
                pgbf1.exponent,
                pgbf1.powers,
                pgbf1.origin,
                pgbf2.exponent,
                (pgbf2.powers[1:ind-1]..., pgbf2.powers[ind]+2, pgbf2.powers[ind+1:end]...),
                pgbf2.origin
            )
            term2 += -0.5 * pgbf2.powers[ind] * (pgbf2.powers[ind]-1) * overlap(
                pgbf1.exponent,
                pgbf1.powers,
                pgbf1.origin,
                pgbf2.exponent,
                (pgbf2.powers[1:ind-1]..., pgbf2.powers[ind]-2, pgbf2.powers[ind+1:end]...),
                pgbf2.origin
            )
        end
        pgbf1.normcoeff * pgbf2.normcoeff * (term0 + term1 + term2)
    end

    kinetic(
            a::ContractedGaussianBasisFunction,
            b::ContractedGaussianBasisFunction
    ) = contract(kinetic, a, b)

    function coulomb(
            pgbf1::PrimitiveGaussianBasisFunction{N},
            pgbf2::PrimitiveGaussianBasisFunction{N},
            pgbf3::PrimitiveGaussianBasisFunction{N},
            pgbf4::PrimitiveGaussianBasisFunction{N}
    ) where N
        r122 = dist2(pgbf1.origin, pgbf2.origin)
        r342 = dist2(pgbf3.origin, pgbf4.origin)

        pcenter = gaussianproductcenter(
            pgbf1.exponent,
            pgbf1.origin,
            pgbf2.exponent,
            pgbf2.origin
        )
        qcenter = gaussianproductcenter(
            pgbf3.exponent,
            pgbf3.origin,
            pgbf4.exponent,
            pgbf4.origin
        )
        rpq2 = dist2(pcenter, qcenter)

        g1 = pgbf1.exponent + pgbf2.exponent
        g2 = pgbf3.exponent + pgbf4.exponent
        delta = .25 * (1/g1 + 1/g2)

        B = ntuple(
            ind->Barray(
                pgbf1.powers[ind],
                pgbf2.powers[ind],
                pgbf3.powers[ind],
                pgbf4.powers[ind],
                pcenter[ind],
                pgbf1.origin[ind],
                pgbf2.origin[ind],
                qcenter[ind],
                pgbf3.origin[ind],
                pgbf4.origin[ind],
                g1,
                g2,
                delta
            ),
            Val(N)
        )

        s = 0.0
        ranges = ntuple(ind->0:pgbf1.powers[ind]+pgbf2.powers[ind]+pgbf3.powers[ind]+pgbf4.powers[ind], Val(N))
        for multi=CartesianIndices(ranges)
            s += prod(B[ind][multi[ind]+1] for ind=1:N) * Fgamma(sum(Tuple(multi)), .25*rpq2/delta)
        end
        return pgbf1.normcoeff*pgbf2.normcoeff*pgbf3.normcoeff*pgbf4.normcoeff*2pi^(2.5)/(g1*g2*sqrt(g1+g2))*exp(-pgbf1.exponent*pgbf2.exponent*r122/g1)*exp(-pgbf3.exponent*pgbf4.exponent*r342/g2)*s
    end

    coulomb(
            a::ContractedGaussianBasisFunction{N},
            b::ContractedGaussianBasisFunction{N},
            c::ContractedGaussianBasisFunction{N},
            d::ContractedGaussianBasisFunction{N},
    ) where N = contract(coulomb, a, b, c, d)

    function Fgamma(m::Int64, x::Float64, SMALL::Float64=1e-12)
        #println("Fgamma($m,$x)")
        x = max(x,SMALL) # Evidently needs underflow protection
        return 0.5 * x^(-m-0.5) * gammainc(m+0.5,x)
    end

    # function gammainc(a::Float64,x::Float64)
    #     # This is the series version of gamma from pyquante. For reasons I don't get, it
    #     # doesn't work around a=1. This works alright, but is only a stopgap solution
    #     # until Julia gets an incomplete gamma function programmed
    #     if abs(a-1) < 1e-3
    #         println("Warning: gammainc_series is known to have problems for a ~ 1")
    #     end
    #     if x < (a+1.0)
    #         #Use the series representation
    #         gam,gln = gser(a,x)
    #     else
    #         #Use continued fractions
    #         gamc,gln = gcf(a,x)
    #         gam = 1-gamc
    #     end
    #     return exp(gln)*gam
    # end
    
    gammainc(a::Float64, x::Float64) = gamma_inc(a, x, 0)[1]::Float64 * gamma(a)

    function gser(a::Float64,x::Float64,ITMAX::Int64=100,EPS::Float64=3e-9)
        # Series representation of Gamma. NumRec sect 6.1.
        gln = (logabsgamma(x))[1]
        if x == 0
            return 0,gln
        end
        ap = a
        delt = s = 1/a
        for i in 1:ITMAX
            ap += 1
            delt *= (x/ap)
            s += delt
            if abs(delt) < abs(s)*EPS
                break
            end
        end
        return s*exp(-x+a*log(x)-gln),gln
    end

    function gcf(a::Float64,x::Float64,ITMAX::Int64=200,EPS::Float64=3e-9,FPMIN::Float64=1e-30)
        #Continued fraction representation of Gamma. NumRec sect 6.1"
        gln = (logabsgamma(x))[1]
        b=x + 1. - a
        c=1. / FPMIN
        d=1. / b
        h=d
        for i in 1:ITMAX
            an=-i*(i-a)
            b=b+2.
            d=an*d+b
            if abs(d) < FPMIN
                d=FPMIN
            end
            c=b+an/c
            if abs(c) < FPMIN
                c=FPMIN
            end
            d= 1. / d
            delt=d*c
            h=h*delt
            if abs(delt-1.) < EPS
                break
            end
        end
        gammcf = exp(-x+a*log(x)-gln)*h
        return gammcf,gln
    end

    function Barray(
            l1::Int64,
            l2::Int64,
            l3::Int64,
            l4::Int64,
            p::Float64,
            a::Float64,
            b::Float64,
            q::Float64,
            c::Float64,
            d::Float64,
            g1::Float64,
            g2::Float64,
            delta::Float64
    )
        Imax = l1 + l2 + l3 + l4 + 1
        B = zeros(Imax)
        for i1=0:(l1+l2)
            for i2=0:(l3+l4)
                for r1=0:div(i1, 2)
                    for r2=0:div(i2, 2)
                        for u=0:(div(i1+i2, 2) - r1 - r2)
                            I = i1 + i2 - 2*(r1 + r2) - u
                            B[I+1] += Bterm(
                                i1,
                                i2,
                                r1,
                                r2,
                                u,
                                l1,
                                l2,
                                l3,
                                l4,
                                p,
                                a,
                                b,
                                q,
                                c,
                                d,
                                g1,
                                g2,
                                delta
                            )
                        end
                    end
                end
            end
        end
        return B
    end


    function Bterm(
            i1::Int64,
            i2::Int64,
            r1::Int64,
            r2::Int64,
            u::Int64,
            l1::Int64,
            l2::Int64,
            l3::Int64,
            l4::Int64,
            Px::Float64,
            Ax::Float64,
            Bx::Float64,
            Qx::Float64,
            Cx::Float64,
            Dx::Float64,
            gamma1::Float64,
            gamma2::Float64,
            delta::Float64
    )
        # THO eq. 2.22
        #print("Bterm($i1,$i2,$r1,$r2,$u,$l1,$l2,$l3,$l4,$Px,$Ax,$Bx,$Qx,$Cx,$Dx,$gamma1,$gamma2,$delta)=")
        val = (-1)^(i2+u) * fB(
                i1,l1,l2,Px,Ax,Bx,r1,gamma1)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2
            ) * fact_ratio2(
                i1 + i2 - 2*(r1+r2),
                u
            ) * (Qx-Px)^(i1+i2-2*(r1+r2)-2*u) / delta^(i1+i2-2*(r1+r2)-u)
        #println("$val")
        return val
    end

    fB(
            i::Int64,
            l1::Int64,
            l2::Int64,
            p::Float64,
            a::Float64,
            b::Float64,
            r::Int64,
            g::Float64
    ) = binomial_prefactor(i, l1, l2, p - a, p - b)*B0(i, r, g)

    B0(
            i::Int64,
            r::Int64,
            g::Float64
    ) = fact_ratio2(i, r) * (4g)^(r-i)

    fact_ratio2(
            a::Int64,
            b::Int64
    ) = factorial(a) * factorial(b) / factorial(a-2b)

    function nuclear_attraction(
            pgbf1::PrimitiveGaussianBasisFunction{N},
            pgbf2::PrimitiveGaussianBasisFunction{N},
            nuccenter::NTuple{N, Float64}
    ) where {N}
        pcenter = gaussianproductcenter(pgbf1.exponent, pgbf1.origin, pgbf2.exponent, pgbf2.origin)
        gamma = pgbf1.exponent + pgbf2.exponent
        rab2 = dist2(pgbf1.origin, pgbf2.origin)
        rcp2 = dist2(nuccenter, pcenter)
        As = Aarrays(
            pgbf1.powers, 
            pgbf2.powers, 
            pcenter .- pgbf1.origin,
            pcenter .- pgbf2.origin,
            pcenter .- nuccenter,
            gamma
        )
        total = 0
        # e.g., over three dimensions, this is the sum over I, J, K, with 
        # I=0:(aI+bI), J=0:(aJ+bJ), K=0:(aK+bK), of
        #   Ax[I+1]*Ay[J+1]*Az[K+1]*Fgamma(I+J+K,rcp2*gamma)
        for multi=CartesianIndices(((0:pgbf1.powers[ind]+pgbf2.powers[ind] for ind=1:N)...,))
            total += (
                prod(As[n][multi[n]+1] for n=1:N)
                * Fgamma(sum(Tuple(multi)), rcp2*gamma)
            )
        end
        val = (
            -2pi * pgbf1.normcoeff * pgbf2.normcoeff
            * exp(-pgbf1.exponent*pgbf2.exponent*rab2/gamma)
            * total / gamma
        )
        return val
    end

    nuclear_attraction(
            cgbf1::ContractedGaussianBasisFunction{N,M},
            cgbf2::ContractedGaussianBasisFunction{N,M},
            nuccenter::NTuple{N, Float64}
    ) where {N, M} = contract(
        (a, b) -> nuclear_attraction(a, b, nuccenter),
        cgbf1, cgbf2
    )
        
    function Aterm(
            i::Int64,
            r::Int64,
            u::Int64,
            l1::Int64,
            l2::Int64,
            ax::Float64,
            bx::Float64,
            cx::Float64,
            gamma::Float64
    )
        term1 = (-1)^(i+u) * binomial_prefactor(i, l1, l2, ax, bx)
        term2 = factorial(i) * cx^(i - 2r - 2u)
        term3 = (1/4/gamma)^(r+u)/factorial(r)/factorial(u)/factorial(i-2r-2u)
        return term1 * term2 * term3
    end
    
    function Aarrays(
            l1s::NTuple{N, Int64},
            l2s::NTuple{N, Int64},
            a::NTuple{N, Float64},
            b::NTuple{N, Float64},
            c::NTuple{N, Float64},
            g::Float64
    ) where {N}
        Imax = ((l1 + l2 + 1 for (l1,l2) = zip(l1s, l2s))...,)
        As = ((zeros(Float64, Imax[n]) for n=1:N)...,)
        for n=1:N
            for i=0:(Imax[n] - 1)
                for r=0:div(i, 2)
                    for u=0:div(i - 2r, 2)
                        I = i - 2r - u + 1
                        As[n][I] += Aterm(
                            i,
                            r,
                            u, 
                            l1s[n], 
                            l2s[n], 
                            a[n], 
                            b[n], 
                            c[n], 
                            g
                        )
                    end
                end
            end
        end
        return As
    end
    

    struct GaussianBasis{M, N} <: AbstractBasis
        cgbfs::NTuple{M, ContractedGaussianBasisFunction{N}}
    end

    Base.length(basis::GaussianBasis{N}) where N = N

    function GaussianBasis(a::Float64, N::Int64)
        cgbfs = ((ContractedGaussianBasisFunction(
            (0,),
            (ind*a/N,),
            (PrimitiveGaussianBasisFunction(
                (N/a)^2,
                (0,),
                (ind*a/N,),
            ),),
            (1.0,))
            for ind=1:N)...,
        )
        GaussianBasis(cgbfs)
    end

    function Laconic.kineticenergyoperator(basis::GaussianBasis{N}) where N
        matrix = zeros(N, N)
        for ind1=1:N
            for ind2=1:ind1
                matrix[ind1, ind2] = kinetic(basis.cgbfs[ind1], basis.cgbfs[ind2])
                matrix[ind2, ind1] = matrix[ind1, ind2]
            end
        end
        Operator("kineticenergy", matrix, basis)
    end

    function nuclearattractionoperator(basis::GaussianBasis{M,N}, nuccenter::Tuple{N, Float64}) where {M,N}
        matrix = zeros(N, N)
        for ind1=1:N
            for ind2=1:N
                matrix[ind1, ind2] = nuclear_attraction(basis.cgbfs[ind1], basis.cgbfs[ind2], nuccenter)
                matrix[ind2, ind1] = matrix[ind1, ind2]
            end
        end
        Operator("nuclearattraction", matrix, basis)
    end

    function coulomboperator(
        basis::CombinedBasis{Tuple{GaussianBasis{M1,N}, GaussianBasis{M2,N}}},
        thresh=1e-7
    ) where {M1, M2, N}
        basis1 = basis.bases[1]
        basis2 = basis.bases[2]
        M = M1*M2
        matrix = zeros(M, M)
        # matrix = SharedArray{Float64}(N, N)
        linind = LinearIndices((1:M1, 1:M2))
        for ind1=1:M1
            # println(Threads.threadid())
            @inbounds for ind2=1:M2
                for ind3=1:M1
                    @simd for ind4=1:M2
                        row = linind[ind2, ind1]
                        col = linind[ind4, ind3]
                        result = coulomb(
                            basis1.cgbfs[ind1],
                            basis1.cgbfs[ind3],
                            basis2.cgbfs[ind2],
                            basis2.cgbfs[ind4]
                        )
                    end
                end
            end
        end
        matrix[matrix .< thresh] .= 0
        Operator("coulomb", sparse(matrix), basis)
    end

    function Laconic.symbolic(basis::GaussianBasis, n::Int64, var::Variable)
        symbolic(basis.cgbfs[n], (var,))
    end

    export PrimitiveGaussianBasisFunction, ContractedGaussianBasisFunction
    export amplitude, overlap, kinetic, coulomb, symbolic
    export nuclear_attraction
    export GaussianBasis
    export coulomboperator, nuclearattractionoperator
end
