module SamplesM
    using Laconic.SystemM
    using Laconic.DisplayM
    using Laconic.Symbolic
    using Laconic.Calculus
    using Laconic
    using Cubature

    using SparseArrays

    function two_electron()
        N1 = 10
        a = 20.
        mass1 = 1
        N2 = 10
        mass2 = 1
        basis1 = DiscreteMomentumBasis(N1, a, mass1)
        pos1 = positionoperator(basis1)
        basis2 = DiscreteMomentumBasis(N2, a, mass2)
        pos2 = positionoperator(basis2)
        basis = kron(basis1, basis2)
        pos = kron(pos1, pos2)
        # TODO: figure out how to represent this as an Operator
        N = N1*N2
        hamiltonian = zeros(N, N)
        for n1=1:N1
            for n2=1:N2
                for n1p=1:N1
                    # Threads.@threads for n2p=1:N2
                    for n2p=1:N2
                        row = N2 * n1 + n2
                        col = N2 * n1p + n2p
                        x1 = Variable("x1")
                        x2 = Variable("x2")
                        psix1 = psix(basis1, n1, x1)
                        psix2 = psix(basis2, n2, x2)
                        psix1p = psix(basis1, n1p, x1)
                        psix2p = psix(basis2, n2p, x2)
                        inversedist = integralidentity(x1, x2)
                        t = Variable("t")
                        # integral = DefiniteIntegral(t, 0, Inf, DefiniteIntegral(x1, 0, a, DefiniteIntegral(x2, 0, a,
                        #    psix1 * psix2 * psix1p * psix2p * 2/sqrt(pi) * Exponential(-t^2*(x1-x2)^2))))
                        # integral = makefinite(integral)
                        w = DefiniteIntegral(x1, 0., a, DefiniteIntegral(x2, 0., a,
                            psix1 * psix2 * psix1p * psix2p * 2/sqrt(pi) * Exponential(-t^2*(x1-x2)^2)
                        ))
                        println(w)
                        t_lin = 2.
                        t_log = 1e4
                        part1 = collapseintegrals(DefiniteIntegral(t, 0., t_lin, w))
                        part2 =  DefiniteIntegral(
                            t,
                            log(t_lin),
                            log(t_log),
                            DefiniteIntegral(w.variable, w.startpoint, w.endpoint,
                                DefiniteIntegral(w.integrand.variable, w.integrand.startpoint, w.integrand.endpoint, w.integrand.integrand * exp(t)))
                        )
                        part2 = collapseintegrals(part2)
                        part3 = pi/t_log^2
                        eval1 = evaluateintegral(part1)
                        println(eval1)
                        eval2 = evaluateintegral(part2)
                        println(eval2)
                        hamiltonian[row,col] = eval1 + eval2 + part3
                    end
                end
            end
        end
    end

    export two_electron
end
