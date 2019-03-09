"""
# module rubidium87

- Julia version:
- Author: david
- Date: 2019-03-02

# Examples

```jldoctest
julia>
```
"""
module rubidium87
    using SparseArrays

    struct Spin
        spin::Rational
        Spin(spin::Rational) = spin.den in (1,2) && spin >= 0 ? new(spin) : error("Invalid spin")
    end
    Spin(spin::Integer) = Spin(Rational(spin))
    convert(::Type{Spin}, x::Integer) = Spin(x)
    convert(::Type{Spin}, x::Rational) = Spin(x)

    function splus(s::Spin)
        Is = collect(1:2*s.spin)
        Js = collect(2:2*s.spin+1)
        Ms = collect(s.spin-1:-1:-s.spin)
        Vs = sqrt.(s.spin*(s.spin+1) .- Ms.*(Ms .+ 1))
        result = sparse(Is, Js, Vs, 2*s.spin+1, 2*s.spin+1)
        result
    end
    sminus(s::Spin) = transpose(splus(s))
    sx(s::Spin) = (splus(s) .+ sminus(s)) ./ 2
    sy(s::Spin) = (splus(s) .- sminus(s)) ./ (2im)
    sz(s::Spin) = sparse(1:2*s.spin+1, 1:2*s.spin+1, s.spin:-1:-s.spin)
    id(s::Spin) = sparse(1:2*s.spin+1, 1:2*s.spin+1, ones(Int(2*s.spin)+1))

    splus(s::Any) = splus(convert(Spin, s))
    sminus(s::Any) = sminus(convert(Spin, s))
    sx(s::Any) = sx(convert(Spin, s))
    sy(s::Any) = sy(convert(Spin, s))
    id(s::Any) = id(convert(Spin, s))

    Ix = kron(sx(3//2), id(1//2), id(0))
    Iy = kron(sy(3//2), id(1//2), id(0))
    Iz = kron(sz(3//2), id(1//2), id(0))
    Sx = kron(id(3//2), sx(1//2), id(0))
    Sy = kron(id(3//2), sy(1//2), id(0))
    Sz = kron(id(3//2), sy(1//2), id(0))
    Lx = kron(id(3//2), id(1//2), id(0))
    Ly = kron(id(3//2), id(1//2), id(0))
    Lz = kron(id(3//2), id(1//2), id(0))
    Jx = Sx .+ Lx
    Jy = Sy .+ Ly
    Jz = Sz .+ Lz
    Fx = Ix .+ Jx
    Fy = Iy .+ Jy
    Fz = Iz .+ Jz

    mu_B = 1.3996245 # Bohr magneton, e hbar/2m_e
    hbar = 0.15915494
    A = 3417.34
    gS = -2.0023193043622
    gL = -0.99999369
    gI = 0.0009951414
    Hhf(Bz::Number) = A * (Ix * Jx + Iy * Jy + Iz * Jz) - mu_B * Bz * (gI*Iz + gS*Sz + gL*Lz)

    export splus, sx, sy, sz, id, Spin, convert
end
