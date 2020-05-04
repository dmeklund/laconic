module SystemM
    using ...Laconic
    using ..Calculus
    using ..Gaussian
    using ..Symbolic
    
    using Combinatorics
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
        x0 = 10.
        sigma = 1.0
        xgrid = Vector(1:cutoff) * a / (cutoff + 1)
        psi0 = normalize(exp.(-((xgrid .- x0) ./ (2*sigma)).^2)) |> Vector{ComplexF64}
        tspan = (0., 30.)
        sol = solve_system(hamiltonian.matrix, basis, psi0, tspan)
        return sol
    end

    function two_particle()
        a = 30.0
        mass1 = mass2 = 1.0
        cutoff = 40
        basis1 = GaussianBasis(a, cutoff)
        basis2 = GaussianBasis(a, cutoff)
        combined_basis = kron(basis1, basis2)
        kinenergy1 = kron(kineticenergyoperator(basis1), identity(basis2))
        kinenergy2 = kron(identity(basis1), kineticenergyoperator(basis2))
        kinenergy = kinenergy1 + kinenergy2
        repulsion = coulomboperator(combined_basis)
        hamiltonian = kinenergy + repulsion
        x1 = 10.
        x2 = 20.
        sigma = 1.0
        xgrid1 = Vector(1:cutoff) * a / (cutoff + 1)
        psi1 = normalize(exp.(-((xgrid1 .- x1) ./ (2*sigma)).^2)) |> Vector{ComplexF64}
        psi2 = normalize(exp.(-((xgrid1 .- x2) ./ (2*sigma)).^2)) |> Vector{ComplexF64}
        psi = kron(psi1, psi2)
        tspan = (0., 30.)
        mat = hamiltonian.matrix
        sol = solve_system(hamiltonian.matrix, combined_basis, psi, tspan)
        return sol
    end

    function apply_at_time(op::AbstractOperator, sln::TimeDependentSolution, t::AbstractFloat)
        state = state_at_time(sln, t)
        apply(op, state)
    end

    flatindex(dims, inds) = LinearIndices(dims)[inds...]

    struct NIdenticalParticleState{T, N}
        statemat::AbstractArray{T, N}
        symmetry::Int8
    end
    getcoeff(
        state::NIdenticalParticleState{T, N},
        bases::NTuple{N, Int64}
    ) where {T, N} = state.statemat[bases...]
    Base.length(state::NIdenticalParticleState{T, N}) where {T, N} = N
    Base.size(state::NIdenticalParticleState{T, N}) where {T, N} = size(state.statemat)
    function asvector(state::NIdenticalParticleState{T, N}) where {T, N}
        dims = size(state)
        statevec = zeros(T, prod(dims))
        # state matrix for identical particles should always be square;
        # i.e., each particle must use the same basis.
        for bases=CartesianIndices(dims)
            coeff = getcoeff(state, Tuple(bases))
            index = flatindex(dims, Tuple(bases))
            statevec[index] = coeff
        end
        statevec
    end
    function setcoeff!(
        state::NIdenticalParticleState{T, N},
        bases::NTuple{N, Int64},
        coeff::T
    ) where {N, T}
        perms = permutations(bases)
        s = collect(1:N)
        next = iterate(perms, s)
        while next !== nothing
            if state.symmetry == -1
                sgn = -2 * parity(s) + 1
            else
                sgn = 1
            end
            permuted, s = next
            state.statemat[permuted...] = sgn * coeff
            next = iterate(perms, s)
        end
    end



    struct NParticleState
        numparticles::Int64
        substates::AbstractVector
    end
    Base.length(state::NParticleState) = state.numparticles
    # Notation is a bit confusing here.
    # If size(a) = (3,3,3), then
    # ( (size(a) for i=1:3)..., ) = ( (3,3,3), (3,3,3), (3,3,3) )
    # which is not what we want - we want to combine into one big tuple of ints,
    # which can be achieved with:
    # ( ((size(a) for i=1:3)...)..., ) = (3,3,3,3,3,3,3,3,3)
    Base.size(state::NParticleState)  = (
        ((size(substate) for substate in state.substates)...)...,
    )
    function getcoeff(
        state::NParticleState,
        bases::NTuple{N, Int64}
    ) where {N}
        # assumes the non-identical particles are uncoorelated
        coeff = 1.0
        offset = 1
        for (substateind, substate) in enumerate(state.substates)
            subbases = bases[offset:offset+length(substate)-1]
            coeff *= getcoeff(substate, subbases)
            offset += length(substate)
        end
        # for (particleind, particlebasis) in enumerate(bases)
        #     state, particleoffset = getsubstate(state, particleind)
        #     coeff *= getcoeff(state, particleoffset)
        # end
        coeff
    end
    function getsubstate(
        state::NParticleState,
        particleind::Int64
    )
        if particleind < 1
            error("Invalid particle index $(particleind)")
        end
        for state in state.substates
            if particleind <= length(state)
                return (state, particleind)
            end
            particleind -= length(state)
        end
        error("Particle index $(particleind) not found in $(state)")
    end
    function asvector(state::NParticleState)
        dims = size(state)
        statevec = zeros(prod(dims))
        for bases=CartesianIndices(dims)
            basestup = convert(Tuple, bases)
            coeff = getcoeff(state, basestup)
            index = flatindex(dims, basestup)
            statevec[index] = coeff
        end
        statevec
    end
    function setcoeff!(
        state::NParticleState,
        groupind::Int64,
        bases::NTuple{N, Int64},
        coeff::T
    ) where {N, T}
        setcoeff!(state.substates[groupind], bases, coeff)
    end

    mutable struct IdenticalParticleBuilder
        type::Type
        N::Int64
        numbasis::Int64
        symmetry::Int8
    end
    Base.length(builder::IdenticalParticleBuilder) = builder.N
    # IdenticalParticleBuilder(
    #         type::Type,
    #         symmetry::Int8) = IdenticalParticleBuilder(0, symmetry)
    function build(builder::IdenticalParticleBuilder)
        NIdenticalParticleState(
            zeros(builder.type, ntuple(n->builder.numbasis, builder.N)),
            builder.symmetry
        )
    end
    

    struct NParticleBuilder
        type::Type
        subbuilders::AbstractVector{IdenticalParticleBuilder}
    end
    NParticleBuilder(type::Type) = NParticleBuilder(type, zeros(type, 0))
    NParticleBuilder() = NParticleBuilder(Float64)
    Base.length(builder::NParticleBuilder) = sum(
        length(subbuilder) for subbuilder in builder.subbuilders
    )
    function addgroup!(
        builder::NParticleBuilder,
        numparticles::Int64,
        numbasis::Int64,
        symmetry::Int64
    )
        push!(
            builder.subbuilders,
            IdenticalParticleBuilder(
                builder.type,
                numparticles,
                numbasis,
                convert(Int8, symmetry)
            )
        )
    end
    function build(builder::NParticleBuilder)
        NParticleState(
            length(builder),
            collect(build(subbuilder) for subbuilder=builder.subbuilders)
        )
    end


    identicalparticlestate(
            type::Type,
            numbasis::Int64, 
            numparticles::Int64,
            symmetry::Int8
    ) = NIdenticalParticleState{type, numparticles}(
        zeros(Float64, ntuple(n->numbasis, numparticles)),
        symmetry
    )

    identicalparticlestate(
        numbasis::Int64,
        numparticles::Int64,
        symmetry::Int8
    ) = identicalparticlestate(Float64, numbasis, symmetry)

    function particlegroup(basis::AbstractBasis, N::Int64, symmetry::Int64)
        numbasis = size(basis)
        dims = ntuple(n->numbasis, N)
        statemat = zeros(ComplexF64, dims)
        for inds=CartesianIndices(dims)

        end
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
    export identicalparticlestate, setcoeff!, getcoeff, asvector
    export NParticleBuilder, NParticleState
    export IdenticalParticleBuilder
    export build, addgroup!, flatindex
end
