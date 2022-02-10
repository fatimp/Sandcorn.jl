porosity(image :: AbstractArray{Bool}) =
    count(iszero, image) / length(image)

function corrfn(image  :: AbstractArray{Bool},
                cutoff :: Integer)
    return Directional.s2(image, false; periodic = true, len = cutoff) |> mean
end

function likelihood_point(x  :: AbstractFloat,
                          xs :: AbstractVector{<:AbstractFloat})
    σ = var(xs)
    μ = mean(xs)

    return -log(sqrt(2π * σ)) - (x - μ)^2/(2σ)
end

function likelihood(cf, size,
                    porosity   :: AbstractFloat,
                    dist       :: Distribution;
                    ncloud     :: Integer = 100)
    cutoff = length(cf)
    cloud = mapreduce(vcat, 1:ncloud) do _
        corrfn(generate_image(size, porosity, dist), cutoff)'
    end

    cloud_points = (cloud[:, n] for n in 1:cutoff)
    res = sum(likelihood_point(x, xs) for (x, xs) in zip(cf, cloud_points))
    return isnan(res) ? -Inf : res
end

function likelihood(image  :: AbstractArray{Bool},
                    dist   :: Distribution;
                    ncloud :: Integer = 100,
                    cutoff :: Integer = typemax(Int))
    likelihood(corrfn(image, cutoff), size(image), porosity(image), dist;
               ncloud = ncloud)
end
