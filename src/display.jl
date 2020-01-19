module DisplayM
    using Interact
    using Plots

    function createSlider(trange, func)
        @manipulate for t=trange
            plot(func(t), ylims=(0,.2))
        end
    end

    function sliderforsoln(soln)
        createSlider(range(soln.tvals[1], soln.tvals[end], length=500), t->abs2.(soln.odesol(t)))
    end

    export createSlider, sliderforsoln
end
