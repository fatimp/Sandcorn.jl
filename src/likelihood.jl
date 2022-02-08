porosity(image :: AbstractArray{Bool}) =
    count(iszero, image) / length(image)

function corrfn(image  :: AbstractArray{Bool},
                cutoff :: Integer)
    s2  = Directional.s2(image, false; periodic = true, len = cutoff) |> mean
    #l2v = Directional.l2(image, false; periodic = true, len = cutoff) |> mean
    #l2s = Directional.l2(image, true;  periodic = true, len = cutoff) |> mean

    return vcat(s2)#, l2v, l2s)
end

function likelihood_point(x  :: AbstractFloat,
                          xs :: AbstractVector{<:AbstractFloat})
    σ = var(xs)
    μ = mean(xs)

    return -log(sqrt(2π * σ)) - (x - μ)^2/(2σ)
end

function likelihood(cf, porosity, image_size, params :: AbstractVector{Tuple{T, T}};
                    ncloud :: Integer = 100,
                    cutoff :: Integer = typemax(Int)) where T <: AbstractFloat
    cloud = mapreduce(vcat, 1:ncloud) do _
        corrfn(generate_image(image_size, porosity, params), cutoff)'
    end

    cloud_points = (cloud[:, n] for n in 1:size(cloud, 2))
    res = sum(likelihood_point(x, xs) for (x, xs) in zip(cf, cloud_points))
    return isnan(res) ? -Inf : res
end

likelihood(image  :: AbstractArray{Bool},
           params :: AbstractVector{Tuple{T, T}};
           ncloud :: Integer = 100,
           cutoff :: Integer = typemax(Int)) where T <: AbstractFloat =
               likelihood(corrfn(image, cutoff), porosity(image), size(image), params;
                          ncloud = ncloud,
                          cutoff = cutoff)
