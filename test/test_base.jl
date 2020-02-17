import Laconic: Basis, State, VectorType, MatrixType, convertToBasis, commutator
import Laconic.SpinM: Spin, sx, sy, sz, id
import LinearAlgebra: dot, I, kron

function exercise2_1()
    increment(theta, psi) = [cos(theta/2), exp(im*psi)*sin(theta/2)]
    decrement(theta, psi) = [-exp(-im*psi)*sin(theta/2), cos(theta/2)]
    basisUpDown = Basis{ComplexF64}("UpDown", 2, ("up", "down"))
    for theta in (0, pi/4, pi)
        for psi in (0, pi/4, pi)
            inc = increment(theta, psi)
            dec = decrement(theta, psi)
            @assert dot(inc, dec) ≈ 0 "$(dot(inc, dec)) !≈ 0"
            @assert dot(inc, inc) ≈ 1 "$(dot(inc, inc)) !≈ 1"
            @assert dot(dec, dec) ≈ 1 "$(dot(dec, dec)) !≈ 1"
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
    @assert commutator(sx_, sy_) == 1.0im*sz_
    @assert (sx_*sx_ + sy_*sy_ + sz_*sz_).matrix ≈ (spin.spin * (spin.spin+1))I
end

function rubidium87()
    nuclear_spin = Spin(3//2)
    electron_spin = Spin(1//2)
    angular_mom  = Spin(0)
    Ix = kron(sx(nuclear_spin), id(electron_spin), id(angular_mom))
    Iy = kron(sy(nuclear_spin), id(electron_spin), id(angular_mom))
    Iz = kron(sz(nuclear_spin), id(electron_spin), id(angular_mom))
    Sx = kron(id(nuclear_spin), sx(electron_spin), id(angular_mom))
    Sy = kron(id(nuclear_spin), sy(electron_spin), id(angular_mom))
    Sz = kron(id(nuclear_spin), sz(electron_spin), id(angular_mom))
    Lx = kron(id(nuclear_spin), id(electron_spin), sx(angular_mom))
    Ly = kron(id(nuclear_spin), id(electron_spin), sy(angular_mom))
    Lz = kron(id(nuclear_spin), id(electron_spin), sz(angular_mom))
    Jx = Sx + Lx
    Jy = Sy + Ly
    Jz = Sz + Lz
    Fx = Ix + Jx
    Fy = Iy + Jy
    Fz = Iz + Jz

    mu_B = 1.3996245 # Bohr magneton, e hbar/2m_e
    hbar = 0.15915494
    A = 3417.34
    gS = -2.0023193043622
    gL = -0.99999369
    gI = 0.0009951414
    Hhf(Bz::Number) = A * (Ix * Jx + Iy * Jy + Iz * Jz) - mu_B * Bz * (gI*Iz + gS*Sz + gL*Lz)
end

rubidium87()
exercise2_1()
exercise3_2()
