module Sandcorn
using Statistics: mean, var
using Distributions
using CorrelationFunctions
using Optim
using Base.Iterators

include("generate.jl")
include("likelihood.jl")
include("optimize.jl")

export generate_image, likelihood, sandcorn_parameters

end # module
