import Laconic: Basis, State, VectorType, MatrixType, convertToBasis, commutator
import Laconic: Spin, sx, sy, sz
using LinearAlgebra: dot

function exercise2_1()
    increment(theta, psi) = [cos(theta/2), exp(im*psi)*sin(theta/2)]
    decrement(theta, psi) = [-exp(-im*psi)*sin(theta/2), cos(theta/2)]
    basisUpDown = Basis{ComplexF64}("UpDown", 2, ("up", "down"))
    for theta in (0, pi/4, pi)
        for psi in (0, pi/4, pi)
            inc = increment(theta, psi)
            dec = decrement(theta, psi)
            @assert dot(inc, dec) ≈ 0
            @assert dot(inc, inc) ≈ 1
            @assert dot(dec, dec) ≈ 1
            xform = hcat(inc, dec)
            basisIncDec = Basis{ComplexF64}(
                "IncDec",
                ("inc", "dec"),
                basisUpDown,
                xform |> MatrixType{ComplexF64}
            )
            state_inc = State{ComplexF64}([1, 0] |> VectorType{ComplexF64}, basisIncDec)
            state_updown = convertToBasis(state_inc, basisUpDown)
            @assert state_updown.vector[1] ≈ inc[1] "$(state_updown.vector[1]) !≈ $(inc[1])"
            @assert state_updown.vector[2] ≈ inc[2] "$(state_updown.vector[2]) !≈ $(inc[2])"
        end
    end
end

function exercise3_2()
    spin = Spin(1//2)
    sx_ = sx(spin)
    sy_ = sy(spin)
    sz_ = sz(spin)
end

exercise2_1()
exercise3_2()
