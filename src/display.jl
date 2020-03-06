module DisplayM
    using Laconic
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
                plot!(p, xrange, func_list[1](t))
            end
            p
        end
    end

    function sliderforsoln(soln::TimeDependentSolution)
        createSlider(
            range(soln.tvals[1], soln.tvals[end], length=500),
            [t->abs2.(soln.odesol(t))]
        )
    end

    function singlebasis(basis, n, t)
        var = Variable("x")
        sum(
            evalexpr(
                symbolic(basis, n, var),
                var,
                t
            ) for n=length(basis)
        )
    end

    function sliderforsoln(soln::TimeDependentSolution{CombinedBasis{T}}, xrange) where T
        var = Variable("x")
        funclist = ((
            t -> abs2.(sum(
                evalexpr(
                    symbolic(basis, n, var),
                    var,
                    xrange
                ) * soln.odesol(t)[n] for n=length(basis)
            )) for basis in soln.basis.bases
        )...,)
        tvals = range(soln.tvals[1], soln.tvals[end], length=500)
        createSlider(tvals, xrange, funclist)
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

    function showbasisfunctions(basis::DiscretePositionBasis, range) where N
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
