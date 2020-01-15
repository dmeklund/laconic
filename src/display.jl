module DisplayM
    using Interact
    using Plots

    function createSlider(trange, func)
        @manipulate for t=trange
            plot(func(t), ylims=(0,.2))
        end
    end

    function sliderforsoln(soln)
        createSlider(1:length(soln.tvals), t->abs2.(soln.coeffs[:,t]))
    end

    export createSlider, sliderforsoln
end
