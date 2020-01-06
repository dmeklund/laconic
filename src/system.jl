module SystemM
    using Laconic
    using Laconic.Calculus
    using DifferentialEquations
    using Laconic.Symbolic
    using SparseArrays
    import FFTW

    hbar = 1.0

    abstract type System end
#     abstract type ElectromagneticSystem <: System end
#     struct TimeIndependentMagneticField <: ElectromagneticSystem
#         SpatialFunction magnetic_field
#     end
#     function hamiltonian(system::TimeIndependentMagneticField)
#     end
    struct SingleParticleSystem
        mass::Number
        basis::AbstractBasis
        cutoff::Integer
        potential::AbstractExpression
    end
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

    struct TimeDependentSolution
        tvals::Vector{Float64}
        coeffs::Matrix{Complex{Float64}}
        basis::AbstractBasis
    end

    function solve_system(hamiltonian, basis, psi0, tspan)
        function func!(psi, p, t)
            1/(im*hbar) * (hamiltonian * psi)
        end
        prob = ODEProblem(func!, psi0, tspan)
        sol = solve(prob)
        TimeDependentSolution(sol.t, hcat(sol.u...), basis)
    end

    function test_solver()
        a = 1.0
        mass = 1.0
        cutoff = 10
        elements = Array(1:cutoff).^2 * pi^2 * hbar^2 / (2 * mass * a^2)
        kineticEnergy = spdiagm(0 => elements)
        xvals = Array(1:cutoff) * a / (cutoff + 1)
        xpos = spdiagm(0 => xvals)  # x operator in the position basis (approx)
        xmom = FFTW.r2r(xpos, FFTW.RODFT00)/(2*(cutoff+1))  # x operator in the momentum basis
        hamiltonian = kineticEnergy + 1000*xmom
        num_basis = 10
        psi0 = zeros(Complex{Float64}, num_basis)
        psi0[1] = 1.0
        tspan = (0., 100.)
        basis = createDiscreteBasis(MomentumBasis(a), num_basis)
        sol = solve_system(hamiltonian, basis, psi0, tspan)
        return sol
    end

    function apply_at_time(op::AbstractOperator, sln::TimeDependentSolution, t::AbstractFloat)
        state = state_at_time(sln, t)
        apply(op, state)
    end

    function state_at_time(sln::TimeDependentSolution, t::AbstractFloat)
        indrange = searchsorted(sln.tvals, t)
        if indrange.start == indrange.stop
            # indicates an exact match was found
            statevec = sln.coeffs[:,indrange.start]
        elseif indrange.stop == 0 || indrange.start > length(sln.t)
            error("$(t) is out of range of the solution")
        else
            coeffs1 = sln.coeffs[indrange.start-1]
            coeffs2 = sln.coeffs[indrange.start]
            tdiff = sln.t[indrange.start] - sln.t[indrange.start-1]
            alpha = (t - sln.t[indrange.start]) / tdiff
            statevec = (1-alpha)*coeffs1 + alpha*coeffs2
        end
        State(statevec, sln.basis)
    end

    export test_solver, apply_at_time
end
