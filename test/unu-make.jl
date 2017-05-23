# Test headers written by `unu make`
using NRRD, FileIO, Colors, AxisArrays, ImageAxes, Unitful, SimpleTraits
using Base.Test

headerpath = joinpath(dirname(@__FILE__), "headers")
writepath = mktempdir()
const mm = u"mm"
const μm = u"μm"
const s = u"s"

@traitimpl TimeAxis{Axis{:t}}

@testset "unu headers" begin

for (file, T, axs, perm) in (("test2duchar.nhdr", UInt8, (5,6), ()),
                             ("test2dushort.nhdr", UInt16, (5,6), ()),
                             ("test2dfloat.nhdr", Float32, (5,6), ()),
                             ("test2drgb.nhdr", RGB{N0f8}, (5,6), ()),
                             ("test2drgb_last.nhdr", RGB{N0f8}, (5,6), (3,1,2)),

                             ("test2drgb_spu_spacing.nhdr",
                              RGB{N0f8},
                              (Axis{:space_1}(NRRD.linspace(0.0,8.0,5)*mm),
                               Axis{:space_2}(NRRD.linspace(0.0,17.5,6)*μm)),
                              ()),

                             ("test2d_spu_spacing_labels.nhdr",
                              UInt8,
                              (Axis{:first}(NRRD.linspace(0.0,8.0,5)*mm),
                               Axis{:second}(NRRD.linspace(0.0,15.0,6)*mm)),
                              ()),

                             ("test2drgb_spu_spacing_labels.nhdr",
                              RGB{N0f8},
                              (Axis{:first}(NRRD.linspace(0.0,8.0,5)*mm),
                               Axis{:second}(NRRD.linspace(0.0,15.0,6)*mm)),
                              ()),

                             ("test3d_units_range_axes.nhdr",
                              UInt8,
                              (Axis{:R}(NRRD.linspace(10.0,14.0,3)*mm),
                               Axis{:A}(NRRD.linspace(12.0,16.0,3)*mm),
                               Axis{:S}(NRRD.linspace(100.0,110.0,3)*mm)),
                              ()),

                             ("test3dlist_origin.nhdr",
                              UInt8,
                              (Axis{:R}(1:3),
                               Axis{:A}(2:4),
                               Axis{:S}(3:5),
                               Axis{:dim_4}(0:1)),
                              ()),

                             ("test3d_time.nhdr",
                              UInt16,
                              (Axis{:space_1}(0:7),
                               Axis{:space_2}(0:9),
                               Axis{:space_3}(0:2),
                               Axis{:time}(0:799)),
                              ()),

                             ("test3d_time_labels.nhdr",
                              UInt16,
                              (Axis{:x}(0:7),
                               Axis{:l}(0:9),
                               Axis{:s}(0:2),
                               Axis{:t}(0:799)),
                              ()),

                             ("test3d_time_labels_units.nhdr",
                              UInt16,
                              (Axis{:x}(NRRD.range(0.0,0.25,8)*μm),
                               Axis{:l}(NRRD.range(0.0,0.25,10)*μm),
                               Axis{:s}(NRRD.range(0.0,2.0,3)*μm),
                               Axis{:time}(0:799)),
                              ()),

                             ("test3d_time_labels_units_s.nhdr",
                              UInt16,
                              (Axis{:x}(NRRD.range(0.0,0.25,8)*μm),
                               Axis{:l}(NRRD.range(0.0,0.25,10)*μm),
                               Axis{:s}(NRRD.range(0.0,2.0,3)*μm),
                               Axis{:time}(NRRD.range(0.0,0.8,800) * s)),
                              ()))
    try
        # Read tests
        filerd = joinpath(headerpath, file)
        Tr, axsr, permr, need_bswap = NRRD.arraytype(filerd)
        @test Tr == T
        @test axsr === axs
        if !(axsr === axs)
            for (axr, ax) in zip(axsr, axs)
                println("axr == ax: ", axr == ax)
                println(axr)
                println(ax)
            end
        end
        @test permr == perm
        if contains(file, "time")
            @test any(istimeaxis, axsr)
        end
        # Write tests
        filewr = joinpath(writepath, file)
        open(query(filerd)) do iord
            skipmagic(iord)
            version, header, keyvals, comments = NRRD.parse_header(iord)
            open(filewr, "w") do iowr
                print(iowr, "NRRD")  # magic bytes
                NRRD.write_header(iowr, version, header, keyvals, comments)
            end
            Tw, axsw, permw, need_bswapw = NRRD.arraytype(filewr)
            @test Tw == Tr
            @test axsw == axsr
            @test permw == permr
        end
    catch err
        warn("failure on file $file")
        rethrow(err)
    end
end

end # @testset "unu headers"
