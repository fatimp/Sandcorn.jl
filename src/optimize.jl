line(x0, x1, x) = x0 + (x1-x0)*x

function target_function(image :: AbstractArray{Bool}, ncloud, cutoff)
    porosity = count(iszero, image) / length(image)

    function target(xs)
        μ, σ = xs
        return -likelihood(image, porosity, μ, σ, ncloud, cutoff)
    end

    return target
end

function sandcorn_parameters(image :: AbstractArray{Bool}; ncloud = 100,
                             μ_bounds = (0.0,  40.0),
                             σ_bounds = (1e-3, 10.0))
    μ_lo, μ_hi = μ_bounds
    σ_lo, σ_hi = σ_bounds

    lower = [μ_lo, σ_lo]
    upper = [μ_hi, σ_hi]

    cutoff = 100

    start = line.(lower, upper, rand(Float64, 2))

    return optimize(target_function(image, ncloud, cutoff), start,
                    ParticleSwarm(lower = lower, upper = upper),
                    Optim.Options(show_trace = true, show_every = 1, iterations = 15))
end
