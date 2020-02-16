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
