using Laconic
using Laconic.Symbolic
using LinearAlgebra

function exercise2_1()
    state_up = Vector([1, 0])
    state_down = Vector([0, 1])
    theta = SSymbol("ϑ")
    phi = SSymbol("ϕ")
    increment(theta, phi) = [
        Cosine(theta/Numeric(2)),
         Exponential(Numeric(1.0im)*phi)*Sine(theta/Numeric(2))
    ]
    decrement(theta, phi) = [
        -Exponential(Numeric(-1.0im)*phi)*Sine(theta/Numeric(2)),
        Cosine(theta/Numeric(2))
    ]
    basisUpDown = Basis{AbstractExpression}("UpDown", 2, ("up", "down"))
    inc = increment(theta, phi)
    dec = decrement(theta, phi)
    println("inc . inc* = $(combineterms(dot(inc, conj(inc))))")
    println("dec . dec* = $(combineterms(dot(dec, conj(dec))))")
    println("inc . dec* = $(combineterms(dot(inc, conj(dec))))")
    @assert combineterms(dot(inc, conj(inc))) == 1
    @assert combineterms(dot(dec, conj(dec))) == 1
    @assert combineterms(dot(inc, conj(dec))) == 0

    xform = hcat(inc, dec)
    @assert is_unitary(xform)
    basisIncDec = Basis{AbstractExpression}(
        "IncDec",
        ("inc", "dec"),
        basisUpDown,
        xform |> MatrixType{AbstractExpression}
    )
end

exercise2_1()
