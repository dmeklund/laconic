# import base: Basis, VectorType, State, convertToBasis
using base

function exercise2_1()
    increment(theta, psi) = [cos(theta/2), exp(im*psi)*sin(theta/2)] |> VectorType
    decrement(theta, psi) = [-exp(-im*psi)*sin(theta/2), cos(theta/2)] |> VectorType
    basisUpDown = Basis("UpDown", 2)
    for theta in (0, pi/4, pi)
        for psi in (0, pi/4, pi)
            inc = increment(theta, psi)
            dec = decrement(theta, psi)
            @assert dot(inc, dec) ≈ 0
            @assert dot(inc, inc) ≈ 1
            @assert dot(dec, dec) ≈ 1
            xform = hcat(inc, dec) |> MatrixType
            basisIncDec = Basis("IncDec", basisUpDown, xform)
            state_inc = State([1, 0] |> VectorType, basisIncDec)
            state_updown = convertToBasis(state_inc, basisUpDown)
            @assert state_updown.vector[1] ≈ inc[1] "$(state_updown.vector[1]) !≈ $(inc[1])"
            @assert state_updown.vector[2] ≈ inc[2] "$(state_updown.vector[2]) !≈ $(inc[2])"
        end
    end
end
