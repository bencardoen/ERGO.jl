# This file is part of ERGO.jl.
#
#  ERGO.jl is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License v3 as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ERGO.jl is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License v3 for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ERGO.jl.  If not, see <https://www.gnu.org/licenses/>.
# Copyright Ben Cardoen, 2018-2021
using Test
using ERGO
using Random
using Statistics
import Images

@testset "ERGO.jl" begin

    @testset "pcor" begin
        Random.seed!(42)
        for i in 1:100
            X = rand(100, 100)
            Y = rand(100, 100)
            p = pearsoncorrelation(X, Y)
            @test p != 1
            p = pearsoncorrelation(X, X)
            @test isapprox(p, 1, atol=0.1)
        end
    end

    @testset "framecorrelation" begin
        # function framecorrelation(images, randomize=false, seed=0, step=1)
        Random.seed!(42)
        for i in 1:100
            frames = [rand(100, 100) for i in 1:100]
            ps = framecorrelation(frames, false, 1, 1)
            qs = framecorrelation(frames, true, 1, 1)
            zs = framecorrelation(frames, true, 4, 1)
            @test qs != zs
            @test ps != qs
            X = ones(100, 100)
            frames = [X .+= (rand(100, 100)/10) for i in 1:100]
            ps = framecorrelation(frames, false, 1, 1)
            @test all(ps .>= 0.85)
            ps = framecorrelation(frames, false, -1, 1)
            @test all(ps .>= 0.5)
        end
    end

    @testset "ifcfiles" begin
        td = mktempdir()
        Random.seed!(42)
        N = 10
        imgs = []
        for i in 1:N
            img = rand(10, 20)
            imgname = joinpath(td, "$i.tif")
            Images.save(imgname, img)
            push!(imgs, imgname)
        end
        #function ifc(images, randomize=false, seed=0, step=1)
        r0 = ifc(imgs, false, 0, 1)
        r1 = ifc(imgs, true, 0, 1)
        r2 = ifc(imgs, false, 0, 2)
        @test r0[1] != r1[1]
        @test r2[1] != r1[1]
    end

    @testset "ccd" begin
        # function check(coord, m, M, b)
#645	    return m+b < coord < M-b
#646	end
        m, M = 1, 10
        b = 0
        for coord in 2:9
            @test check(coord, m, M, b)
        end
        b = 1
        for coord in 3:8
            @test check(coord, m, M, b)
        end
    end

    @testset "norm" begin
        xs = [1 2 3]
        nm = normimg(xs)
        @test nm == [1/3 2/3 3/3]
        nm = normimg(xs, true)
        @test nm == [0 1/2 4/4]
        Random.seed!(32)
        for _ in 1:100
            img = rand(10, 10)
            ni = normimg(img, false)
            @test all(0 .<= abs.(ni) .<= 1)
            img = 0.5 .- rand(10, 10)
            cni = normimg(img, true)
            @test all(0 .<= cni .<= 1)
        end
        m = normimg(zeros(10,10))
        @test iszero(m)
        m = normimg(zeros(10,10,10))
        @test iszero(m)
        m = normimg(rand(10,10,10))
        @test !isnothing(m)
    end

    @testset "gmpr" begin
        x, y = 0, 0.0
        r = geom_pair(x, y)
        @test r == 0
        x = 1.0
        r = geom_pair(x, y)
        @test r == 0
        y = 2.0
        r = geom_pair(x, y)
        @test isapprox(r, √2.0)
    end

    @testset "hmpr" begin
        x, y = zeros(10), zeros(10)
        res = harmonic_pair(x, y)
        @test all(res .== 0)
    end

    @testset "mcc" begin
        fp = 1
        tp = 6
        fn = 2
        tn = 3
        mc = mcc(tp, tn, fp, fn)
        @test isapprox(mc, 0.47809144373375745)
    end

    @testset "cc" begin
        vertices = [[0,0],[0,1],[1,0], [2,2], [2,1]]
        adjmat = buildAdjMat(vertices)
        graph3 = buildGraph(adjmat, 0.8)
        graph2 = buildGraph(adjmat, 1.1)
        graph1 = buildGraph(adjmat, 3)
        N = size(vertices)[1]
        cc3 = connectedComponents(N, graph3)
        cc2 = connectedComponents(N, graph2)
        cc1 = connectedComponents(N, graph1)
        @test size(cc3)[1] == N
        @test size(cc2)[1] == 2
        @test size(cc1)[1] == 1
    end

    @testset "gm" begin
        Random.seed!(42)
        b = 1 .+ rand(100)
        _g = gm(b)
        @test _g <= Statistics.mean(b)
    end

    @testset "peaks" begin
        Random.seed!(42)
        a = rand(100, 100)
        p = peaks(a)
        @test !isnothing(p)
    end

    @testset "gmsm" begin
        Random.seed!(42)
        b = 1 .+ rand(100)
        am = Statistics.mean(b)
        gm, gs = gmsm(b)
        @test gm <= am
        @test isapprox(gm ,exp(Statistics.mean(log.(b))), atol=0.01)
        @test isapprox(gs, exp(Statistics.std(log.(b))), atol=0.01)
        a = [0 1 2]
        b = copy(a)
        g, s = gmsm(a)
        @test g!=0
        @test all(a .== b)
    end

    @testset "tomask" begin
        Random.seed!(42)
        X = 10
        Y = 20
        Z = 2
        img1 = rand(X,Y,Z)
        i2 = copy(img1)
        msk = tomask(img1)
        @test all(i2 .== img1)
        @test all(msk[i2 .>0].== 1)
        @test all(msk[i2 .== 0].== 0)
    end

    @testset "aszero" begin
        Random.seed!(42)
        X = 10
        Y = 20
        Z = 2
        img1 = rand(X,Y,Z)
        img2 = copy(img1)
        img0 = aszero(img1)
        @test all(img0 .== zero(eltype(img1)))
        @test all(img1 .== img2)
    end
end
