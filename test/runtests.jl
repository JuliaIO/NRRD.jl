using FileIO, ColorTypes, Unitful, AxisArrays, ImageAxes
using Base.Test
# using Images, ColorTypes, FixedPointNumbers, FileIO

include("unu-make.jl")

@testset "NRRD" begin
    workdir = joinpath(tempdir(), "Images")
    writedir = joinpath(workdir, "write")
    if !isdir(workdir)
        mkdir(workdir)
    end
    if !isdir(writedir)
        mkdir(writedir)
    end

    @testset "Gray, raw" begin
        println("small"); flush(STDOUT)
        img = load(joinpath(dirname(@__FILE__), "io", "small.nrrd"))
        flush(STDOUT)
        @test eltype(img) == Float32
        @test ndims(img) == 3
        @test size(img) == (10,10,10)
        outname = joinpath(writedir, "small.nrrd")
        # save(outname, img)
        # imgc = load(outname)
        # @test img == imgc
    end

    @testset "Units" begin
        println("units"); flush(STDOUT)
        img = load(joinpath(dirname(@__FILE__), "io", "units.nhdr"))
        axsv = axisvalues(img)
        @test step(axsv[1]) == 0.1u"mm"
        @test step(axsv[2]) == 0.2u"mm"
        @test step(axsv[3]) == 1u"mm"
        axsn = axisnames(img)
        @test axsn == (:L, :P, :S)
    end

    @testset "Compressed (gzip)" begin
        println("compressed"); flush(STDOUT)
        img = load(joinpath(dirname(@__FILE__), "io", "smallgz.nrrd"))
        @test ndims(img) == 3
        @test eltype(img) == Float32
        outname = joinpath(writedir, "smallgz.nrrd")
        # save(outname, img)
        # imgc = load(outname)
        # @test img == imgc
    end

    @testset "Time is 4th dimension" begin
        println("small_time"); flush(STDOUT)
        img = load(joinpath(dirname(@__FILE__), "io", "small_time.nrrd"))
        @test timedim(img) == 4
        outname = joinpath(writedir, "small_time.nrrd")
        # save(outname, img)
        # imgc = load(outname)
        # @test img == imgc
    end
end
