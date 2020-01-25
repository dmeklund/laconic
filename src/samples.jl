module SamplesM
    using Laconic.SystemM
    using Laconic.DisplayM
    using Laconic.Symbolic
    using Laconic

    function two_electron()
        N1 = 10
        a = 20
        mass1 = 1
        N2 = 10
        mass2 = 1
        basis1 = DiscreteMomentumBasis(N1, a, mass1)
        basis2 = DiscreteMomentumBasis(N2, a, mass2)
        basis = kron(basis1, basis2)
    end

    export two_electron
end
