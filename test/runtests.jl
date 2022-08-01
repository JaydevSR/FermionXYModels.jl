using FermionXYModels
using LinearAlgebra
using Test

@testset "FermionBasis" begin
    @test length(FermionBasis(3)) == 8
    @test collect(FermionBasis(4)) == vec(collect(Iterators.product([-1, 1], [-1, 1], [-1, 1], [-1, 1])))
    @test_throws ArgumentError FermionBasis(4; states=(-1, 0, 1))
    @test_throws ArgumentError FermionBasis(0)
    @test_throws ArgumentError FermionBasis(-1)
end

@testset "Models" begin
    L = 20
    m1 = FermionXYChain(; L=L, J=1.0, h=1.0, gamma=1.0, start=:vacuum)
    m2 = FermionXYChain(; L=L, J=1, h=1, gamma=1, start=:vacuum)

    fnames = fieldnames(FermionXYChain)
    @test all([getfield(m1, f) for f in fnames] .== [getfield(m2, f) for f in fnames])
    @test_throws ArgumentError FermionXYChain(; L=0, J=1.0, h=1.0, gamma=1.0, start=:rand)
    @test_throws ArgumentError FermionXYChain(; L=L, J=0, h=1.0, gamma=1.0, start=:rand)
    @test_throws ArgumentError FermionXYChain(; L=-20, J=1.0, h=1.0, gamma=1.0, start=:rand)
    @test_throws ArgumentError FermionXYChain(; L=L, J=1.0, h=0, gamma=-2.0, start=:rand)

    m3 = FermionIsingChain(; L=L, J=1.0, h=1.0, start=:vacuum)
    m4 = FermionXYChain(; L=L, J=1.0, h=1.0, gamma=1.0, start=:vacuum)
    @test all([getfield(m3, f) for f in fnames] .== [getfield(m4, f) for f in fnames])

    m4 = FermionXXChain(; L=L, J=1.0, h=1.0, start=:vacuum)
    m5 = FermionXYChain(; L=L, J=1.0, h=1.0, gamma=0, start=:vacuum)
    @test all([getfield(m4, f) for f in fnames] .== [getfield(m5, f) for f in fnames])
end

@testset "Correlations and Probabilities" begin
    p = 0
    L = 10
    for sites in FermionBasis(L)
        p += det(probability_matrix(sites; L=L, J=1.0, h=1.0, gamma=1.0))
    end
    @test isapprox(p, 1.0)
    
    model = FermionXYChain(;L=L, J=1.0, h=1.2, gamma=0.5)
    G_mat = correlation_matrix(model)
    p_mat = probability_matrix(model)

    @test G_mat == correlation_matrix(;L=L, J=1.0, h=1.2, gamma=0.5)
    @test p_mat == probability_matrix(model.sites; L=L, J=1.0, h=1.2, gamma=0.5)
end

