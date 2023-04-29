using MoM_Basics
using Test

@testset "MoM_Basics.jl" begin

    proj_path = joinpath(@__DIR__, "..")
    @testset ".nas formats" begin

        filenames = readdir(joinpath(proj_path, "meshfiles"))
        for fn in filenames
            @show fn
            filename = joinpath(joinpath(proj_path, "meshfiles", fn))
            meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
            @test true

            show_memory_time()
            @test true
        end

    end

    @testset "Triangle, RWG" begin
        filename = joinpath(proj_path, "meshfiles/Tri.nas")
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true
        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData)
        @test true
        show_memory_time()
        @test true
    end

    @testset "Terahedron, PWC and SWG" begin

        filename = joinpath(proj_path, "meshfiles/Tetra.nas")
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :PWC)
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :SWG)
        @test true

        setGeosPermittivity!(geosInfo, 2(1. - 0.001im))
        @test true

        show_memory_time()
        @test true

    end

    
    @testset "Hexadron, PWC and RBF" begin

        filename = joinpath(proj_path, "meshfiles/Hexa.nas")
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :PWC)
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :RBF)
        @test true

        setGeosPermittivity!(geosInfo, 2(1. - 0.001im))
        @test true

        show_memory_time()
        @test true

    end


    @testset "Tetra + Hexadron, PWC" begin

        filename = joinpath(proj_path, "meshfiles/TetraHexa.nas")
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :PWC)
        @test true

        setGeosPermittivity!(geosInfo, 2(1. - 0.001im))
        @test true

        show_memory_time()
        @test true

    end

    @testset "Tri + Tetra, RWG + SWG, RWG + PWC" begin

        filename = joinpath(proj_path, "meshfiles/TriTetra.nas")
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData; sbfT = :RWG, vbfT = :SWG)
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData; sbfT = :RWG, vbfT = :PWC)
        @test true

        setGeosPermittivity!(geosInfo, 2(1. - 0.001im))
        @test true

        show_memory_time()
        @test true

    end

    @testset "Tri + Hexa, RWG + RBF, RWG + PWC" begin

        filename = joinpath(proj_path, "meshfiles/TriHexa.nas")
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData; sbfT = :RWG, vbfT = :RBF)
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData; sbfT = :RWG, vbfT = :PWC)
        @test true

        setGeosPermittivity!(geosInfo, 2(1. - 0.001im))
        @test true

        show_memory_time()
        @test true

    end

    @testset "Tri + Tetra + Hexa, RWG + PWC" begin

        filename = joinpath(proj_path, "meshfiles/TriTetraHexa.nas")
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true
        
        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData; sbfT = :RWG, vbfT = :PWC)
        @test true

        setGeosPermittivity!(geosInfo, 2(1. - 0.001im))
        @test true

        show_memory_time()
        @test true

    end

    rm("results"; force = true, recursive = true)

end
