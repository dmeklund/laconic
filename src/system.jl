module SystemM
    using Laconic
    using Laconic.Calculus
    using Laconic.Gaussian
    using Laconic.Symbolic
    using LinearAlgebra
    using DifferentialEquations
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

    function propagator(hamiltonian::Operator, delta_t::Number)
        Operator("propagator", exp(-im*delta_t/hbar * hamiltonian.matrix), hamiltonian.basis)
    end

    function propagator_magnus(hamiltonian, num_terms)
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

    struct TimeDependentSolution{T <: AbstractBasis}
        tvals::Vector{Float64}
        coeffs::Matrix{Complex{Float64}}
        basis::T
        odesol
    end

    function Laconic.symbolic(soln::TimeDependentSolution, t::Float64, var::Variable)
        coeffs = soln.odesol(t)
        sum(coeffs[n] * symbolic(soln.basis, n, var) for n=1:length(soln.basis))
    end

    function Laconic.symbolic(soln::TimeDependentSolution{CombinedBasis{T}}, basisind::Int, t::Float64, var::Variable) where {T}
        coeffs = soln.odesol(t)
        sum(coeffs[m] * symbolic(soln.basis, m, basisind, var) for m=1:length(soln.basis))
    end

    function solve_system(hamiltonian, basis, psi0, tspan)
        function func!(dpsi, psi, p, t)
            dpsi[:] = 1/(im*hbar) * (hamiltonian * psi)
        end
        prob = ODEProblem(func!, psi0, tspan)
        sol = solve(prob)
        TimeDependentSolution(sol.t, hcat(sol.u...), basis, sol)
    end

    function test_solver()
        a = 30.0
        mass = 1.0
        cutoff = 100
        # elements = Array(1:cutoff).^2 * pi^2 * hbar^2 / (2 * mass * a^2)
        # basis = DiscretePositionBasis(cutoff, a, mass)
        basis = GaussianBasis(a, cutoff)
        kineticEnergy = kineticenergyoperator(basis)
        # xop = positionoperator(basis) # x operator in the momentum basis
        hamiltonian = kineticEnergy #+ xop
        # hamiltonian = xop
        # psi0 = zeros(Complex{Float64}, cutoff)
        # psi0[:] = eigen(hamiltonian.matrix).vectors[:,1] + (rand(80).-.5)* .1
        # psi0 /= sum(psi0 .* psi0)
        # println(psi0)
        # psi0[10] = 1.0
        x0 = 15
        sigma = 1.0
        xgrid = Vector(1:cutoff) * a / (cutoff + 1)
        psi0 = normalize(exp.(-((xgrid .- x0) ./ (2*sigma)).^2)) |> Vector{ComplexF64}
        tspan = (0., 100.)
        sol = solve_system(hamiltonian.matrix, basis, psi0, tspan)
        return sol
    end

    function two_particle()
        a = 30.0
        mass1 = mass2 = 1.0
        cutoff = 20
        basis1 = GaussianBasis(a, cutoff)
        basis2 = GaussianBasis(a, cutoff)
        combined_basis = kron(basis1, basis2)
        kinenergy1 = kron(kineticenergyoperator(basis1), identity(basis2))
        kinenergy2 = kron(identity(basis1), kineticenergyoperator(basis2))
        kinenergy = kinenergy1 + kinenergy2
        repulsion = coulomboperator(combined_basis)
        hamiltonian = repulsion # .1*kinenergy + repulsion
        x1 = 14.
        x2 = 16.
        sigma = 1.0
        xgrid = Vector(1:cutoff) * a / (cutoff + 1)
        psi1 = normalize(exp.(-((xgrid .- x1) ./ (2*sigma)).^2)) |> Vector{ComplexF64}
        psi2 = normalize(exp.(-((xgrid .- x2) ./ (2*sigma)).^2)) |> Vector{ComplexF64}
        psi = kron(psi1, psi2)
        tspan = (0., 10.)
        sol = solve_system(-hamiltonian.matrix, combined_basis, psi, tspan)
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

    export test_solver, apply_at_time, propagator
    export two_particle
    export TimeDependentSolution
end
