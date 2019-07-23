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

    function is_orthonormal(matrix::MatrixType{T}) where T
        product = transpose(matrix) * matrix
        product == I
    end

    include("basis.jl")
    include("operator.jl")
    include("spin.jl")
    include("Symbolic.jl")
end
