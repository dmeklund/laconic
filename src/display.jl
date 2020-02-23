module DisplayM
    using Laconic
    using Laconic.Gaussian
    using Laconic.Calculus
    using Laconic.Symbolic
    using Interact
    using Plots

    function createSlider(trange, func)
        @manipulate for t=trange
            plot(func(t), ylims=(0,.2))
        end
    end

    function sliderforsoln(soln)
        createSlider(
            range(soln.tvals[1], soln.tvals[end], length=500),
            t->abs2.(soln.odesol(t))
        )
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
end
