using WrightOmega
using Test

@testset "WrightOmega.jl" begin
    z_vals = [
        # R1
        4im
        -1 + 3.5im
        1 + 3.5im
        # R2
        -4im
        -1 - 3.5im
        1 - 3.5im
        # R3
        -5 + 0im
        -5
        -3 + 2im
        -3 + 0.9999π * im
        -3 + π * im # uses regularization
        -3 - 0.999π * im # uses regularization
        # R4
        0 + 0im
        0
        0.5im
        -1 + 0im
        -1
        1 + im
        1 - im
        1 + π * im
        1 - π * im
        # R5
        -5 + 3im / 2
        -2 + (π + 3 / 4)im
        # R6
        -5 - 3im / 2
        -2 - (π + 3 / 4)im
        # R7
        5 + 0im
        5
        5 + 5im
        ]
    conv = [Float32, identity, big]
    for z in z_vals, T in conv
        zT = T(real(z)) + T(imag(z)) * im
        w = wrightomega(zT)
        @test w + log(w) ≈ zT
    end
    @test wrightomega(NaN) === NaN
    @test wrightomega(NaN + NaN * im) === NaN + NaN * im
end
