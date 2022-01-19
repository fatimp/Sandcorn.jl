function generate_image(size,
                        porosity :: AbstractFloat,
                        μ        :: AbstractFloat,
                        σ        :: AbstractFloat)
    radius_distribution = Truncated(Normal(μ, σ), 0, Inf)
    image  = zeros(Bool, size)
    cur_pores = prod(size)
    pores     = porosity * cur_pores

    indices    = CartesianIndices(image)
    fidx, lidx = first(indices), last(indices)
    uidx       = oneunit(fidx)

    while true
        radius  = rand(radius_distribution)
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
