# to create new sys_image:
# import PackageCompiler
# symbols = [:Plots, :WebIO, :QuadGK, :DifferentialEquations, :Cubature, :Dierckx, :Reactive, :Mux, :GR, :Interact]
# PackageCompiler.create_sysimage(symbols; sysimage_path="sys_laconic.so", precompile_statements_file="precompile.jl")
# Then start julia with
# julia --sysimage sys_laconic.so


using Plots
using WebIO
using QuadGK
using DifferentialEquations
using Cubature
using Dierckx
using Reactive
using Mux
using GR
using Interact

p = Plots.plot(rand(5), rand(5))
display(p)
webio_serve(page("/", req -> node(:p, "Hello")))