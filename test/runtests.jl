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
        img = load(joinpath(dirname(@__FILE__), "io", "small.nrrd"))
        flush(STDOUT)
        @test eltype(img) == Float32
        @test ndims(img) == 3
        @test size(img) == (10,10,10)
        outname = joinpath(writedir, "small.nrrd")
        isfile(outname) && rm(outname)
        save(outname, img)
        imgc = load(outname)
        @test img == imgc
        @test axisnames(imgc) == axisnames(img)
        @test axisvalues(imgc) == axisvalues(img)
    end

    @testset "Units" begin
        img = load(joinpath(dirname(@__FILE__), "io", "units.nhdr"))
        axsv = axisvalues(img)
        @test step(axsv[1]) == 0.1u"mm"
        @test step(axsv[2]) == 0.2u"mm"
        @test step(axsv[3]) == 1u"mm"
        axsn = axisnames(img)
        @test axsn == (:L, :P, :S)
        outname = joinpath(writedir, "units.nrrd")
        isfile(outname) && rm(outname)
        save(outname, img)
        imgc = load(outname)
        @test img == imgc
        axsv = axisvalues(imgc)
        @test step(axsv[1]) == 0.1u"mm"
        @test step(axsv[2]) == 0.2u"mm"
        @test step(axsv[3]) == 1u"mm"
        axsn = axisnames(imgc)
        @test axsn == (:L, :P, :S)
    end

    @testset "Compressed (gzip)" begin
        img = load(joinpath(dirname(@__FILE__), "io", "smallgz.nrrd"))
        @test ndims(img) == 3
        @test eltype(img) == Float32
        outname = joinpath(writedir, "smallgz.nrrd")
        isfile(outname) && rm(outname)
        save(outname, img)
        imgc = load(outname)
        @test img == imgc
    end

    @testset "Time is 4th dimension" begin
        img = load(joinpath(dirname(@__FILE__), "io", "small_time.nrrd"))
        @test timedim(img) == 4
        outname = joinpath(writedir, "small_time.nrrd")
        isfile(outname) && rm(outname)
        save(outname, img)
        imgc = load(outname)
        @test img == imgc
        @test axisnames(img)[4] == axisnames(imgc)[4] == :time
        @test axisvalues(img) == axisvalues(imgc)
    end
end

nothing
