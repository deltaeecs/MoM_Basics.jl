inputBasicParameters(;frequency = rand()*1e9)

# setPrecision
@test all([Float32, Float64]) do T
    setPrecision!(T)
    Precision.FT == T
end

# setParams
@test isapprox(Params.λ_0, C_0/Params.frequency)
@test isapprox(Params.K_0, 2π/Params.λ_0)