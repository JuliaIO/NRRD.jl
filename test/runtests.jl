using FileIO, FixedPointNumbers, ColorTypes, Unitful, AxisArrays, ImageAxes
using Base.Test

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

    @testset "eltype" begin
        for (elty, name) in ((UInt8, "uint8"),
                             (UInt8, "u8"),
                             (Gray{N0f8}, "gray_u8"),
                             (RGB{N0f8}, "rgb_u8"),
                             (Gray{N0f16}, "gray_u16"),
                             (Gray{N4f12}, "gray_u12"),
                             (Gray{N4f12}, "gray_u12_hex"))
            fn = joinpath(dirname(@__FILE__), "io", "eltype_$name.nrrd")
            img = load(fn)
            @test eltype(img) == elty
            @test size(img) == (3,5)
        end
    end

    gc()  # to close any mmapped files
    rm(workdir, recursive=true)
end

include("readremote.jl")

nothing
