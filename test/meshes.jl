proj_path = joinpath(@__DIR__, "..")
@testset ".nas formats" begin

    filenames = readdir(joinpath(proj_path, "meshfiles"))
    for fn in filenames
        filename = joinpath(joinpath(proj_path, "meshfiles", fn))
        meshData, εᵣs   =  getMeshData(filename; meshUnit=:mm);
        @test true

        show_memory_time()
        @test true
    end

end