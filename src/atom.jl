module Atom
    using Laconic.Gaussian

    function sto3g(atomno::Integer, electronno::Integer)
        if atomno == 1 && electronno == 1
            # 1s orbital
            powers = (0, 0, 0)
            origin = (0., 0., 0.)
            ContractedGaussianBasisFunction(
                powers,
                origin,
                (
                    PrimitiveGaussianBasisFunction(3.4252509099999999, powers, origin),
                    PrimitiveGaussianBasisFunction(0.62391373000000006, powers, origin),
                    PrimitiveGaussianBasisFunction(0.16885539999999999, powers, origin)
                ),
                (0.15432897000000001, 0.53532813999999995, 0.44463454000000002)
            )
        end
    end

    export sto3g
end
