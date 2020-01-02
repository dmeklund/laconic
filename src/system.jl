module SystemM
    using Laconic.Calculus
    using DifferentialEquations
    abstract type System end
#     abstract type ElectromagneticSystem <: System end
#     struct TimeIndependentMagneticField <: ElectromagneticSystem
#         SpatialFunction magnetic_field
#     end
#     function hamiltonian(system::TimeIndependentMagneticField)
#     end
    function propagator(hamiltonian, num_terms)
        result = 0
        for ind in 1:num_terms
            result += propagator_element(hamiltonian, ind)
        end
    end

    function propagator_element(hamiltonian, index, t0, t)
        if index == 1
            t1 = Variable("t1")
            return -im/hbar * DefiniteIntegral(t1, t0, t, hamiltonian(t1))
        elseif index == 2
            t1 = Variable("t1")
            t2 = Variable("t2")
            return -1/(2*hbar^2) * DefiniteIntegral(t1, t0, t, DefiniteIntegral(t2, t0, t1, commutator(hamiltonian(t1), hamiltonian(t2))))
        elseif index == 3
            t1 = Variable("t1")
            t2 = Variable("t2")
            t3 = Variable("t3")
            return im/(6*hbar^3) * DefiniteIntegral(t1, t0, t,
                DefiniteIntegral(t2, t0, t1,
                    DefiniteIntegral(t3, t0, t2, (commutator(hamiltonian(t1), commutator(hamiltonian(t2), hamiltonian(t3))))))
            )
        else
            error("Only up to three Magnus expansion elements supported")
        end
    end

    function solve_system(hamiltonian, basis, psi0, tspan)
        function func!(dpsi, psi, p, t)
            dpsi = 1/(im*hbar) * (hamiltonian * psi)
        end
        prob = ODEProblem(func!, psi0, tspan)
        sol = solve(prob)
    end
end
