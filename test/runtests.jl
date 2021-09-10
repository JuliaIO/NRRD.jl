using FileIO, ImageCore, Unitful, AxisArrays, ImageAxes, ImageMetadata
using Test, Base.CoreLogging

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
        flush(stdout)
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
        # Check that FixedPoint types get properly encoded
        outname = joinpath(writedir, "small14.nrrd")
        img = rand(Gray{N2f14}, 5, 5, 2)
        save(outname, img)
        imgr = load(outname)
        @test eltype(imgr) == Gray{N2f14}
        @test imgr == img
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
        img = AxisArray(rand(N0f8, 6, 5, 3, 2), :x, :y, :z, :time)
        save(outname, img)
        imgr = load(outname)
        @test axisnames(imgr)[4] == :time
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
            outname = joinpath(writedir, "elt.nrrd")
            save(outname, img)
            imgr = load(outname)
            @test eltype(imgr) == elty
        end
    end

    @testset "bool" begin
        # Bool is not one of the supported eltypes in the `type` field of
        # http://teem.sourceforge.net/nrrd/format.html.
        # Make sure we don't do something wrong
        fn = joinpath(writedir, "eltype_bool.nrrd")
        img = rand(Bool, 5, 6)
        save(fn, img)
        imgr = load(fn)
        @test imgr == img
        # We don't insist on eltype(imgr) == Bool
    end

    @testset "Mmapped" begin
        fn = joinpath(writedir, "volume.nrrd")
        save(fn, zeros(UInt8,500,500,500))
        v = load(fn)
        @test size(v) == (500,500,500)
        @test all(x->x==0, v)
    end

    @testset "Endian mmapped" begin
        fn = joinpath(dirname(@__FILE__), "io", "bswap.nrrd")
        img = load(fn; mmap=true)
        @test img == [255 254; 1 2]
    end

    @testset "Fiji compatibility" begin
        for (name, ps) in (("fiji-stack-pixels-16bit.nrrd", (1,1,1)),
                           ("fiji-stack-microns-16bit.nrrd", (0.4μm, 0.4μm, 2.0μm)))
            fn = joinpath(dirname(@__FILE__), "io", name)
            v = load(fn)
            @test size(v) == (50,40,30)
            @test pixelspacing(v) == ps
        end
    end

    @testset "Header only" begin
        img = load(joinpath(dirname(@__FILE__), "io", "small.nrrd"))
        outname = joinpath(writedir, "small.nhdr")
        isfile(outname) && rm(outname)
        props = Dict("datafile"=>joinpath(writedir,"small.raw"))
        save(outname, img; props=props)
        imgc = load(outname)
        @test img == imgc
        @test axisnames(imgc) == axisnames(img)
        @test axisvalues(imgc) == axisvalues(img)
    end

    @testset "N0f16 ImageMeta" begin
        img = ImageMeta(rand(N0f16, 5, 7), info=false)
        outname = joinpath(writedir, "imagemeta.nhdr")
        save(outname, img)
        imgc = load(outname)
        @test img == imgc

        # With units
        img = AxisArray(img, (:y, :x), (1.0u"mm", 1.2u"mm"))
        save(outname, img)
        imgc = load(outname)
        @test img == imgc
        @test pixelspacing(imgc) == (1.0u"mm", 1.2u"mm")
        img = ImageMeta(img, info=false)
        save(outname, img)
        imgc = load(outname)
        @test img == imgc
        @test AxisArrays.HasAxes(img) isa AxisArrays.HasAxes{true}
        @test pixelspacing(imgc) == (1.0u"mm", 1.2u"mm")
    end

    @testset "Size warning" begin
        for (T, Traw) in ((N0f8, UInt8), (N0f16, UInt16))
            img = rand(T, 5, 7, 3)
            outname = joinpath(writedir, "imagetrunc.nhdr")
            props = Dict("datafile"=>joinpath(writedir, "imagetrunc.raw"))
            save(outname, img, props=props)
            # Create a too-small data array
            open(props["datafile"], "w") do io
                write(io, rand(Traw, 5, 7, 2))
            end
            # Test that it can be read but produces a warning
            logger = Test.TestLogger()
            with_logger(logger) do
                imgr = load(outname)
                @test size(imgr) == (5, 7, 2)
            end
            record = logger.logs[1]
            @test occursin("header indicates an array size (5, 7, 3), but the file size is consistent with at most (5, 7, 2)", record.message)
        end
    end

    @testset "Extended eltype" begin
        img = load(joinpath(dirname(@__FILE__), "io", "test_helix.nrrd"))
        @test axisnames(img)[end-2:end] == (:R, :A, :S)
    end

    GC.gc()  # to close any mmapped files
    try
        rm(workdir, recursive=true)
    catch
    end
end

include("readremote.jl")

nothing
