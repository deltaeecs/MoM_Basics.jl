using MoM_Basics
using Test

@testset "MoM_Basics.jl" begin

    @testset "Params" begin
        include("params.jl")
    end

    @testset "Meshes" begin
        include("meshes.jl")
    end

    @testset "Basis functions" begin
        include("basis_functions.jl")
    end
    rm("results"; force = true, recursive = true)

end
