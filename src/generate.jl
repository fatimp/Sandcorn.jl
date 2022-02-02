function rand_non_negative(distribution)
    r = rand(distribution)
    return r < 0 ? rand_non_negative(distribution) : r
end

"""
    generate_image(size, porosity, normal_param)

Generate an image of overlapping `n`-dimensional balls with radii
sampled from a truncated mixture of normal distributions with
parameters `normal_param` in the form of `[(μ₀, σ₀), (μ₁, σ₁),
…]`. The image is generated with specified `size` and `porosity`.
"""
function generate_image(size,
                        porosity     :: AbstractFloat,
                        normal_param :: AbstractVector{Tuple{T, T}}) where T <: AbstractFloat
    radius_distribution = MixtureModel(Normal, normal_param)
    # Truncated seems to be broken,
    # see https://github.com/JuliaStats/Distributions.jl/issues/1501
    #radius_distribution = Truncated(radius_distribution, 0, Inf)

    image     = zeros(Bool, size) |> BitArray
    cur_pores = prod(size)
    pores     = porosity * cur_pores

    indices    = CartesianIndices(image)
    fidx, lidx = first(indices), last(indices)
    uidx       = oneunit(fidx)

    while true
        radius  = rand_non_negative(radius_distribution)
        iradius = radius |> ceil |> Int
        center  = rand(indices)

        from = max(center - iradius*uidx, fidx)
        to   = min(center + iradius*uidx, lidx)

        old_pores = count(iszero, image[from:to])

        for idx in from:to
            r = sqrt(sum(x^2 for x in Tuple(idx - center)))
            if r <= radius
                image[idx] = 1
            end
        end

        new_pores = count(iszero, image[from:to])
        cur_pores += new_pores - old_pores

        if cur_pores < pores
            return image
        end
    end
end
