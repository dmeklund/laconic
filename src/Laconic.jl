# using Revise
# push!(LOAD_PATH, pwd())

"""
# module base.jl

- Julia version: 1.0
- Author: deklund
- Date: 2019-06-09

# Examples

```jldoctest
julia>
```
"""
module Laconic
    import LinearAlgebra: transpose, I
    MatrixType{T} = AbstractArray{T, 2}
    VectorType{T} = AbstractArray{T, 1}

    include("basis.jl")
    include("Symbolic.jl")
    include("operator.jl")
    include("calculus.jl")
    include("system.jl")
    include("display.jl")
    include("spin.jl")
    include("samples.jl")

    export MatrixType, VectorType, is_unitary, commutator
end
