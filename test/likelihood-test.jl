@testset "Test likelihood function" begin
    param = [(5.0, 1.0)]
    image = generate_image((1000, 1000), 0.15, param)
    # 1.18 is a factor of safety since all this thing is purely random
    # and can never achieve its "true" maximum value.
    max_likelihood = 1.18*Sandcorn.likelihood(image, param; cutoff = 20, ncloud = 150)

    for _ in 1:5
        new_param = [(5 + 1 - 2rand(Float64), 1 + 0.5 - rand(Float64))]
        shifted_likelihood = Sandcorn.likelihood(image, new_param; cutoff = 20, ncloud = 150)
        @test shifted_likelihood < max_likelihood
    end
end
