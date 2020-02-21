using Laconic.Gaussian

function test_pgbf()
    s = PrimitiveGaussianBasisFunction(1.0)
    px = PrimitiveGaussianBasisFunction(1.0, (1,0,0), (0.,0.,0.))
    amp1 = amplitude(s, (0.,0.,0.))
    amp2 = amplitude(px, (0.,0.,0.))
    # print("Amplitude1 ($(amp1)) should be .71270547")
    # print("Amplitude2 ($(amp2)) should be 0")
    @assert isapprox(amplitude(s, (0.,0.,0.)), 0.71270547)
    @assert isapprox(amplitude(px, (0.,0.,0.)), 0)
end
test_pgbf()

function test_cgbf()
    c = ContractedGaussianBasisFunction((0,0,0), (0.0, 0.0, 0.0), (1.0,1.0))
    # push!(c,1,1)
    @assert isapprox(amplitude(c, (0.,0.,0.)), 0.71270547)
    c2 = ContractedGaussianBasisFunction(
        (0,0,0),
        (0., 0., 0.),
        (1.0, 0.2),
        (0.5, .2)
    )
    @assert isapprox(overlap(c2,c2),1)
end
test_cgbf()

function test_kinetic()
    s = PrimitiveGaussianBasisFunction(1.0)
    c = ContractedGaussianBasisFunction((0,0,0), (0.0,0.0,0.0), (1.0, 1.0))
    @assert isapprox(amplitude(c, (0.,0.,0.)), 0.71270547)
    # @assert isapprox(kinetic(1.,0.,0.,0.,0,0,0,1.,0.,0.,0.,0,0,0),2.9530518648229536)
    @assert isapprox(kinetic(s,s),1.5) "$(kinetic(s,s)) !≈ 1.5"
    @assert isapprox(kinetic(c,c),1.5) "$(kinetic(c,c)) !≈ 1.5"
end
test_kinetic()

using Laconic.Gaussian: fB, B0, fact_ratio2, Bterm
function test_two_terms()
    @assert fB(0,0,0,0.0,0.0,0.0,0,2.0) == 1
    @assert fB(0,0,0,1.0,1.0,1.0,0,2.0) == 1
    @assert fB(0,0,0,0.0,0.0,0.0,0,2.0 ) == 1
    @assert fB(1,0,1,0.0,0.0,0.0,0,2.0 ) == 0.125
    @assert B0(0,0,2.0) == 1
    @assert fact_ratio2(0,0) == 1
    @assert Bterm(0,0,0,0,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,0.25)==1
    @assert Bterm(0,1,0,0,0,0,0,0,1,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,0.25)==0
end
test_two_terms()

function test_coul1()
    s = PrimitiveGaussianBasisFunction(1.0)
    px = PrimitiveGaussianBasisFunction(1.0, (1,0,0), (0.,0.,0.))
    @assert coulomb(s,s,s,px)==0 # 0
    @assert isapprox(coulomb(s,s,px,px), 0.9403159725793305 )
end
test_coul1()
