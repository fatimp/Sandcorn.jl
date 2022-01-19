corrfn(image :: AbstractArray{Bool}) = Directional.s2(image, false; periodic = true) |> mean

function likelihood_point(s  :: AbstractFloat,
                          cs :: AbstractVector{T}) where T <: AbstractFloat
    σ = var(cs)
    μ = mean(cs)

    return log(1/(sqrt(2π * σ))) - (s - μ)^2/(2σ)
end

function likelihood(image :: AbstractArray{Bool}, porosity, μ, σ, ncloud, cutoff = typemax(Int))
    cloud = reduce(hcat, (generate_image(size(image), porosity, μ, σ) |> corrfn for _ in 1:ncloud))
    cf    = corrfn(image)

    len = min(cutoff, length(cf))
    res = sum(likelihood_point(cf[n], cloud[n, :]) for n in 1:len)
    return isnan(res) ? -Inf : res
end
