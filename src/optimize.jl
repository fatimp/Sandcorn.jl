line(x0, x1, x) = x0 + (x1-x0)*x

function array2params(xs :: AbstractArray{<:AbstractFloat})
    len = length(xs)
    @assert len |> iseven

    μs = xs[begin:len÷2]
    σs = xs[len÷2 + 1:end]
    return zip(μs, σs) |> collect
end

function target_function(image  :: AbstractArray{Bool},
                         ncloud :: Integer,
                         cutoff :: Integer)
    cf = corrfn(image, cutoff)
    p  = porosity(image)

    return xs -> -likelihood(cf, p, size(image), array2params(xs);
                             ncloud = ncloud,
                             cutoff = cutoff)
end

"""
~~~~
sandcorn_parameters(image; nmodes     = 1,
                           ncloud     = 100,
                           cutoff     = 100,
                           iterations = 15,
                           μ_bounds   = (0.0,  40.0),
                           σ_bounds   = (1e-3, 10.0))
~~~~

Find a distribution of radii in an image which consists of partially
overlapping grains (like an slice of a sandstone). `nmodes` is an
assumed number of modes in the distribution. The resulting
distribution can be constructed as
`Truncated(MixtureModel(Normal, sandcorn_parameters(…)), 0, Inf)`.
"""
function sandcorn_parameters(image      :: AbstractArray{Bool};
                             nmodes     :: Integer = 1,
                             ncloud     :: Integer = 100,
                             cutoff     :: Integer = 100,
                             iterations :: Integer = 15,
                             μ_bounds = (0.0,  40.0),
                             σ_bounds = (1e-3, 10.0))
    μ_lo, μ_hi = μ_bounds
    σ_lo, σ_hi = σ_bounds

    μ_lower = fill(μ_lo, nmodes)
    σ_lower = fill(σ_lo, nmodes)
    μ_upper = fill(μ_hi, nmodes)
    σ_upper = fill(σ_hi, nmodes)

    lower = vcat(μ_lower, σ_lower)
    upper = vcat(μ_upper, σ_upper)
    start = line.(lower, upper, rand(Float64, 2nmodes))

    opt = optimize(target_function(image, ncloud, cutoff), start,
                   ParticleSwarm(lower = lower, upper = upper),
                   Optim.Options(show_trace = true,
                                 show_every = 1,
                                 iterations = iterations))
    return opt |> Optim.minimizer |> array2params
end
