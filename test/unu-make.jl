# Test headers written by `unu make`
using NRRD, FileIO, Colors, AxisArrays, ImageAxes, Ranges, Unitful, SimpleTraits
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
                             ("test2drgb.nhdr", RGB{U8}, (5,6), ()),
                             ("test2drgb_last.nhdr", RGB{U8}, (5,6), (3,1,2)),

                             ("test2drgb_spu_spacing.nhdr",
                              RGB{U8},
                              (Axis{:space_1}(Ranges.linspace(0mm,8mm,5)),
                               Axis{:space_2}(Ranges.linspace(0μm,17.5μm,6))),
                              ()),

                             ("test2d_spu_spacing_labels.nhdr",
                              UInt8,
                              (Axis{:first}(Ranges.linspace(0mm,8mm,5)),
                               Axis{:second}(Ranges.linspace(0mm,15mm,6))),
                              ()),

                             ("test2drgb_spu_spacing_labels.nhdr",
                              RGB{U8},
                              (Axis{:first}(Ranges.linspace(0mm,8mm,5)),
                               Axis{:second}(Ranges.linspace(0mm,15mm,6))),
                              ()),

                             ("test3d_units_range_axes.nhdr",
                              UInt8,
                              (Axis{:R}(Ranges.linspace(10mm,14mm,3)),
                               Axis{:A}(Ranges.linspace(12mm,16mm,3)),
                               Axis{:S}(Ranges.linspace(100mm,110mm,3))),
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
                              (Axis{:x}(Ranges.range(0μm,0.25μm,8)),
                               Axis{:l}(Ranges.range(0μm,0.25μm,10)),
                               Axis{:s}(Ranges.range(0μm,2μm,3)),
                               Axis{:time}(0:799)),
                              ()),

                             ("test3d_time_labels_units_s.nhdr",
                              UInt16,
                              (Axis{:x}(Ranges.range(0μm,0.25μm,8)),
                               Axis{:l}(Ranges.range(0μm,0.25μm,10)),
                               Axis{:s}(Ranges.range(0μm,2μm,3)),
                               Axis{:time}(Ranges.range(0s,0.8s,800))),
                              ()))
    try
        # Read tests
        filerd = joinpath(headerpath, file)
        Tr, axsr, permr, need_bswap = NRRD.arraytype(filerd)
        @test Tr == T
        @test axsr === axs
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
