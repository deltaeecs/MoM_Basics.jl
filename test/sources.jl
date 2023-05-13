
@testset "PlaneWave" begin

    setPrecision!(Float64)

    plw1 = PlaneWave{Float64}(π, π, 0.)
    plw2 = PlaneWave(π, π, 0.)
    @test plw1 == plw2

    e = sourceEfield(plw1, [0, 0, 0])
    @test norm(e) == 1

    h = sourceHfield(plw1, [0, 0, 0])
    @test norm(h) == 1/η_0

end


@testset "MagneticDipole" begin

    md = MagneticDipole()

    # set phase
    p0 = md.Iml
    update_phase!(md, π)
    p1 = md.Iml
    @test p0 == -p1

    # add phase
    add_phase!(md, π)
    p2 = md.Iml
    @test isapprox(p0, p2)

    # update orient
    update_orient!(md, [0, π/2, 0])
    @test md.orient == [0, π/2, 0]

    # sourceEfield
    elc = sourceLocalEfield(md, [ 1., 0., 0.])
    update_orient!(md, [π/2, 0, 0])
    egl = sourceEfield(md, [0., 1., 0.])
    @test isapprox(md.l2gRot*elc, egl)

    egls = sourceEfield([md, md], [0., 1., 0.])
    @test isapprox(egls, 2egl)

    # sourceFarEfield
    efarlc = sourceLocalFarEfield(md, r̂θϕInfo(π/2, π/2))
    phase_zxz = [2π*rand(), π*rand(), 0]
    update_orient!(md, [phase_zxz[1], phase_zxz[2], phase_zxz[3]])
    efargl = sourceFarEfield(md, r̂θϕInfo(π/2-phase_zxz[2], π/2+phase_zxz[1]))
    @test isapprox(efarlc, efargl)

    # update_orient!(md, [0, 0, 0])
    efargls = sourceFarEfield([md, md], r̂θϕInfo(π/2-phase_zxz[2], π/2+phase_zxz[1]))
    @test isapprox(efargls, 2efargl)

    # radiationDirectionCoeff
    @test radiationDirectionCoeff(md, θϕInfo(0, 0)) ≈ 0
    @test radiationDirectionCoeff(md, θϕInfo(π/2, π/2)) ≈ 1.5

end

@testset "AntennaArray" begin

    ary = antennaArray((2, 2), [0., 0., 0.]; sourceConstructer = MagneticDipole, 
                        sourceT = MagneticDipole, sourceorientlc = [π/2, π/2, 0.], coefftype = :taylor)
    setdiffArray!(ary, 1)

    e1      =   sourceLocalEfield(ary, [1., 1., 1.])
    efar1   =   sourceLocalFarEfield(ary, r̂θϕInfo([1., 1., 0.]))
    update_orient!(ary; aryorient = [0., π/2, 0.], sourceorientlc = [π/2, π/2, 0.])
    e2      =   sourceEfield(ary, ary.l2gRot*[1., 1., 1.])
    efar2   =   sourceFarEfield(ary, r̂θϕInfo(ary.l2gRot*[1., 1., 0.]))

    l2gmat  =   eulerRotationMat(ary.orient...)

    @test isapprox(l2gmat*e1, e2)

    @test isapprox(norm(efar1), norm(efar2))
    
end