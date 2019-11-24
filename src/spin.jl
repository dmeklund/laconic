module SpinM
    import LinearAlgebra: transpose
    import SparseArrays: sparse
    import Laconic: Basis, MatrixType, Operator

    struct Spin
        spin::Rational
        function Spin(spin::Rational)
            if spin.den in (1,2) && spin >= 0
                new(spin)
            else
                error("Invalid spin")
            end
        end
        Spin(spin::Integer) = Spin(spin//1)
    end
    dicke(s::Spin) = Basis{ComplexF64}("Dicke", (2*s.spin).num+1)
    function splus(s::Spin)
        Is = 1:2*s.spin
        Js = 2:2*s.spin+1
        Ms = s.spin-1:-1:-s.spin
        Vs = sqrt.(s.spin*(s.spin+1) .- Ms.*(Ms.+1)) .+ 0.0im
        basis = dicke(s)
        Operator{ComplexF64}("splus", sparse(Is, Js, Vs, 2*s.spin+1, 2*s.spin+1), basis)
    end
    sminus(s::Spin) = transpose(splus(s), "sminus")
    sx(s::Spin) = Operator{ComplexF64}("sx", (splus(s) + sminus(s)) / 2.0)
    sy(s::Spin) = Operator{ComplexF64}("sy", (splus(s) - sminus(s)) / 2.0im)
    sz(s::Spin) = Operator{ComplexF64}(
        "sz",
        sparse(1:2*s.spin+1, 1:2*s.spin+1, s.spin:-1:-s.spin), # |> VectorType{ComplexF64}),
        dicke(s)
    )
    id(s::Spin) = Operator{ComplexF64}(
        "id",
        sparse(1:2*s.spin+1, 1:2*s.spin+1, ones((2*s.spin).num+1) |> VectorType{ComplexF64}),
        dicke(s)
    )

    # Pauli matrices
    σ_x = [0 1; 1 0]
    σ_y = [0 -im; im 0]
    σ_z = [1 0; 0 -1]

    export Spin, sx, sy, sz, splus, sminus, id
    export σ_x, σ_y, σ_z
end
