module DisplayM
    using ...Laconic
    using Laconic.Gaussian
    using Laconic.Calculus
    using Laconic.Symbolic
    using Laconic.SystemM
    using Interact
    using Plots

    function createSlider(trange, xrange, func_list)
        @manipulate for t=trange
            p = plot()
            for func in func_list
                plot!(p, xrange, func(t))
            end
            p
        end
    end

    function singlebasis(basis::AbstractBasis, n, t)
        var = Variable("x")
        sum(
            evalexpr(
                symbolic(basis, n, var),
                var,
                t
            ) for n=1:length(basis)
        )
    end

    function sliderforsoln(soln::TimeDependentSolution, xrange)
        var = Variable("x")
        basis = soln.basis
        func = t -> abs.(evalexpr(
            symbolic(soln, t, var),
            var,
            xrange
        ))
        tvals = range(soln.tvals[1], soln.tvals[end], length=500)
        createSlider(tvals, xrange, (func,))
    end

    function sliderforsoln(soln::TimeDependentSolution{CombinedBasis{T}}, xrange) where T
        var = Variable("x")
        funclist = ((
            t -> abs2.(evalexpr(
                symbolic(soln, basisind, t, var),
                var,
                xrange
            )) for basisind=1:length(soln.basis.bases)
        )...,)
        tvals = range(soln.tvals[1], soln.tvals[end], length=500)
        createSlider(tvals, xrange, funclist)
    end

    function plotattime(soln::TimeDependentSolution{CombinedBasis{T}}, time, xrange) where T
        var = Variable("x")
        p = plot()
        for basisind=1:length(soln.basis.bases)
            plot!(
                p, 
                xrange, 
                abs2.(evalexpr(
                    symbolic(soln, basisind, time, var),
                    var,
                    xrange
                ))
            )
        end
    end

    function showbasisfunctions(basis::GaussianBasis{N}, range) where N
        x = Variable("x")
        p = plot()
        for ind=1:N
            sym = symbolic(basis.cgbfs[ind], (x,))
            func = convertToFunction(sym, x)
            plot!(p, range, func.(range))
        end
        p
    end

    function showbasisfunctions(basis::DiscretePositionBasis, range)
        x = Variable("x")
        p = plot()
        for ind=1:basis.N
            sym = symbolic(basis, ind, x)
            func = convertToFunction(sym, x)
            plot!(p, range, func.(range))
        end
        p
    end

    export createSlider, sliderforsoln, showbasisfunctions
    export singlebasis
end
