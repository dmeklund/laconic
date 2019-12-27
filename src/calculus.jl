#=
calculus:
- Julia version: 1.2.0
- Author: david
- Date: 2019-12-26
=#
module Calculus
    struct IntegrationVariable
        label::String
    end
    struct DefiniteIntegral{T1, T2, T3}
        variable::IntegrationVariable
        startpoint::T1
        endpoint::T2
        integrand::T3
    end
    export IntegrationVariable, DefiniteIntegral
end