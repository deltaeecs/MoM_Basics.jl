using MoM_Basics
using Test

@testset "MoM_Basics.jl" begin

    @testset ".Nas Formats" begin

        filenames = [   "sphere_r1m_metal_300MHz.nas", "radomeTetra15mm.nas", 
                        "radomeHexa15mm.nas", "missileHT15mm.nas",
                        "sphereShellTetra50mm.nas", "sphereShellHexa50mm.nas",
                        "sphereShellHT50mm.nas"]
        for fn in filenames
            @show fn
            filename = joinpath("../meshfiles", fn)
            meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
            @test true

            show_memory_time()
            @test true
        end

    end

    @testset "Triangle, RWG" begin
        filename = "../meshfiles/sphere_r1m_metal_300MHz.nas"
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true
        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData)
        @test true
        show_memory_time()
        @test true
    end

    @testset "Terahedron, PWC and SWG" begin

        filename = "../meshfiles/radomeTetra15mm.nas"
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :PWC)
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :SWG)
        @test true

        setGeosPermittivity!(geosInfo, 1 + 1im)
        @test true

        show_memory_time()
        @test true

    end

    
    @testset "Hexadron, PWC and RBF" begin

        filename = "../meshfiles/radomeHexa15mm.nas"
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :PWC)
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :RBF)
        @test true

        setGeosPermittivity!(geosInfo, 1 + 1im)
        @test true

        show_memory_time()
        @test true

    end


    @testset "Tetra + Hexadron, PWC" begin

        filename = "../meshfiles/missileHT15mm.nas"
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData, vbfT = :PWC)
        @test true

        setGeosPermittivity!(geosInfo, 1 + 1im)
        @test true

        show_memory_time()
        @test true

    end

    rm("results"; force = true, recursive = true)

end
