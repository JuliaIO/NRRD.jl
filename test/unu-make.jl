# Test headers written by `unu make`
using NRRD, ImageAxes, Unitful

headerpath = joinpath(dirname(@__FILE__), "headers")
const mm = u"mm"
const μm = u"μm"

for (file, T, axs, perm) in (("test2duchar.nhdr", UInt8, (5,6), ()),
                             ("test2dushort.nhdr", UInt16, (5,6), ()),
                             ("test2dfloat.nhdr", Float32, (5,6), ()),
                             ("test2drgb.nhdr", RGB{U8}, (5,6), ()),
                             ("test2drgb_last.nhdr", RGB{U8}, (5,6), (3,1,2)),

                             ("test2drgb_spu_spacing.nhdr",
                              RGB{U8},
                              (Axis{:row}(0mm:2mm:8mm),
                               Axis{:col}(0.0μm:3.5μm:17.50001μm)),
                              ()),

                             ("test2d_spu_spacing_labels.nhdr",
                              UInt8,
                              (Axis{:first}(0mm:2mm:8mm),
                               Axis{:second}(0mm:3mm:15mm)),
                              ()),

                             ("test2drgb_spu_spacing_labels.nhdr",
                              RGB{U8},
                              (5,6),
                              (Axis{:first}(0mm:2mm:8mm),
                               Axis{:second}(0mm:3mm:15mm)),
                              ()),

                             ("test3d_units_range_axes.nhdr",
                              UInt8,
                              (Axes{:R}(10mm:2mm:14mm),
                               Axes{:A}(12mm:2mm:18mm),
                               Axes{:S}(100mm:5mm:110mm)),
                              ()),

                             ("test3dlist_origin.nhdr",
                              UInt8,
                              (Axes{:R}(1:4),
                               Axes{:A}(2:5),
                               Axes{:S}(3:6),
                               Axis{:dim_4}(1:2)),
                              ()),

                             ("test3d_time.nhdr",
                              UInt16,
                              (Axes{:row}(1:8),
                               Axes{:col}(1:10),
                               Axes{:page}(1:3),
                               Axis{:time}(1:800)),
                              ()),

                             ("test3d_time_labels.nhdr",
                              UInt16,
                              (Axes{:x}(1:8),
                               Axes{:l}(1:10),
                               Axes{:s}(1:3),
                               Axis{:t}(1:800)),
                              ()),

                             ("test3d_time_labels_units.nhdr",
                              UInt16,
                              (Axes{:x}(0μm:0.25μm:1.76μm),
                               Axes{:l}(0μm:0.25μm:2.26μm),
                               Axes{:s}(0μm:2μm:4μm),
                               Axis{:t}(1:800)),
                              ()),

                             ("test3d_time_labels_units_s.nhdr",
                              UInt16,
                              (Axes{:x}(0μm:0.25μm:1.76μm),
                               Axes{:l}(0μm:0.25μm:2.26μm),
                               Axes{:s}(0μm:2μm:4μm),
                               Axis{:t}(0s:0.8s:639.21s)),
                              ()))
    Tr, axsr, permr, need_bswap = arraytype(joinpath(headerpath, file))
    @test Tr == T
    @test axsr == axs
    @test permr == perm
    if contains(file, "time")
        @test timeaxis(axsr) != nothing
    end
end
