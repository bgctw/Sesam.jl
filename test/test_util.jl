test_smoothstep = (x_step, dx, a, b) -> begin
    # at the center between edges
    c = (a+b)/2
    @test smoothstep(x_step,x_step,dx,a,b) == c
    # before and at edge0 return a
    @test smoothstep(x_step - dx -0.0001,x_step,dx,a,b) == a
    @test smoothstep(x_step - dx,x_step,dx,a,b) == a
    # after and at edge1 return b
    @test smoothstep(x_step + dx +0.0001,x_step,dx,a,b) == b
    @test smoothstep(x_step + dx, x_step,dx,a,b) == b
    if a <= b
        # before x_step return between a and center
        @test a <= smoothstep(x_step - dx/2,x_step,dx,a,b) <= c
        # after x_step return between center and b
        @test c <= smoothstep(x_step + dx/2,x_step,dx,a,b) <= b
    else
        # before x_step return between a and center
        @test a >= smoothstep(x_step - dx/2,x_step,dx,a,b) >= c
        # after x_step return between center and b
        @test c >= smoothstep(x_step + dx/2,x_step,dx,a,b) >= b
    end
    # monotonous
    xs = range(x_step-dx,x_step+dx,length=10)
    ys = smoothstep.(xs,x_step,dx,a,b)
    dys = [ys[i+1] - ys[i] for i in 1:(length(ys)-1)]
    @test all(sign.(dys) .== sign(b-a))
    true
end

@testset "smoothstep" begin
    x = 0; x_step = 0; dx = 1; a = 0; b = 1
    test_smoothstep(0,1,0,1)
    test_smoothstep(0,1,2,4)
    test_smoothstep(0,1,4,2)
    test_smoothstep(0,1,4,4)
    #
    test_smoothstep(3,1,0,1)
    test_smoothstep(3,1,2,4)
    test_smoothstep(3,1,4,2)
    test_smoothstep(3,1,4,4)
    #
    test_smoothstep(-3,1,0,1)
    test_smoothstep(-3,1,2,4)
    test_smoothstep(-3,1,4,2)
    test_smoothstep(-3,1,4,4)
    #
    @test_throws AssertionError smoothstep(3,0,0,0,1)
end;

i_tmp = () -> begin
    # using CairoMakie
    ts = -0.2:0.01:2
    scatter(ts, smoothstep.(ts,0.5,0.2) .* (1 .- smoothstep.(ts, 1, 0.2)))
end
