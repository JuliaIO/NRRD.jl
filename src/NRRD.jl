module NRRD

using Mmap, Printf
# Packages needed to return the possible range of element types
using FixedPointNumbers, Colors, ColorVectorSpace, StaticArrays, Quaternions
# Other packages
using AxisArrays, ImageAxes, Unitful, MappedArrays
using FileIO
import Libz
import FixedPointNumbers

using Colors: AbstractGray
using AxisArrays: HasAxes

const string2type = Dict(
    "signed char" => Int8,
    "int8" => Int8,
    "int8_t" => Int8,
    "uchar" => UInt8,
    "unsigned char" => UInt8,
    "uint8" => UInt8,
    "uint8_t" => UInt8,
    "short" => Int16,
    "short int" => Int16,
    "signed short" => Int16,
    "signed short int" => Int16,
    "int16" => Int16,
    "int16_t" => Int16,
    "ushort" => UInt16,
    "unsigned short" => UInt16,
    "unsigned short int" => UInt16,
    "uint16" => UInt16,
    "uint16_t" => UInt16,
    "int" => Int32,
    "signed int" => Int32,
    "int32" => Int32,
    "int32_t" => Int32,
    "uint" => UInt32,
    "unsigned int" => UInt32,
    "uint32" => UInt32,
    "uint32_t" => UInt32,
    "longlong" => Int64,
    "long long" => Int64,
    "long long int" => Int64,
    "signed long long" => Int64,
    "signed long long int" => Int64,
    "int64" => Int64,
    "int64_t" => Int64,
    "ulonglong" => UInt64,
    "unsigned long long" => UInt64,
    "unsigned long long int" => UInt64,
    "uint64" => UInt64,
    "uint64_t" => UInt64,
    "float16" => Float16,
    "float" => Float32,
    "double" => Float64
)

# The opposite of string2type
type2string(::Type{Float16}) = "float16"
type2string(::Type{Float32}) = "float"
type2string(::Type{Float64}) = "double"
type2string(::Type{T}) where {T<:Integer} = lowercase(string(T.name.name))
type2string(::Type{Normed{T,f}}) where {T<:Unsigned,f} = type2string(T)
type2string(::Type{T}) where {T} = type2string(eltype(T), T)
type2string(::Type{T}, ::Type{T}) where {T} = error("type $T unrecognized")
type2string(::Type{T1}, ::Type{T2}) where {T1,T2} = type2string(T1)

const space2axes = Dict(
    "right-anterior-superior" => (3,(:R,:A,:S)),
    "ras" => (3,(:R,:A,:S)),
    "left-anterior-superior" => (3,(:L,:A,:S)),
    "las" => (3,(:L,:A,:S)),
    "left-posterior-superior" => (3,(:L,:P,:S)),
    "lps" => (3,(:L,:P,:S)),
    "right-anterior-superior-time" => (4,(:R,:A,:S,:time)),
    "rast" => (4,(:R,:A,:S,:time)),
    "left-anterior-superior-time" => (4,(:L,:A,:S,:time)),
    "last" => (4,(:L,:A,:S,:time)),
    "left-posterior-superior-time" => (4,(:L,:P,:S,:time)),
    "lpst" => (4,(:L,:P,:S,:time)),
    "scanner-xyz" => (3,(:scannerx,:scannery,:scannerz)),
    "scanner-xyz-time" => (4,(:scannerx,:scannery,:scannerz,:time)),
    "3d-right-handed" => (3,(:xrcs,:yrcs,:zrcs)),
    "3d-left-handed" => (3,(:xlcs,:ylcs,:zlcs)),
    "3d-right-handed-time" => (4,(:xrcs,:yrcs,:zrcs,:time)),
    "3d-left-handed-time" => (4,(:xlcs,:ylcs,:zlcs,:time)),
)

const axes2space = Dict(
    (:R,:A,:S) => "right-anterior-superior",
    (:L,:A,:S) => "left-anterior-superior",
    (:L,:P,:S) => "left-posterior-superior",
    (:R,:A,:S,:time) => "right-anterior-superior-time",
    (:L,:A,:S,:time) => "left-anterior-superior-time",
    (:L,:P,:S,:time) => "left-posterior-superior-time",
    (:scannerx,:scannery,:scannerz) => "scanner-xyz",
    (:scannerx,:scannery,:scannerz,:time) => "scanner-xyz-time",
    (:xrcs,:yrcs,:zrcs) => "3d-right-handed",
    (:xlcs,:ylcs,:zlcs) => "3d-left-handed",
    (:xrcs,:yrcs,:zrcs,:time) => "3d-right-handed-time",
    (:xlcs,:ylcs,:zlcs,:time) => "3d-left-handed-time")

# We put these in a dict so that we don't eval untrusted
# strings. Please submit PRs to add to this list if you need
# additional unit support.
const unit_string_dict = Dict("" => 1, "m" => u"m", "mm" => u"mm", "s" => u"s",
                              "um" => u"μm", "μm" => u"μm", "microns" => u"μm",
                              "pixel" => 1)

struct QString end                 # string with quotes around it: "mm"
VTuple{T} = Tuple{Vararg{T}}  # space-delimited tuple: 80 150
struct PTuple{T} end               # parenthesis-delimited tuple: (80,150)
struct StringPTuple{T} end         # string or PTuple{T}

struct IntFloat end

# This should list anything that DOESN'T parse to a string
const parse_type = Dict(
    # basic
    "dimension"=>Int,
    "block size"=>Int,
    "blocksize"=>Int,
    "min"=>Float64,
    "max"=>Float64,
    "old min"=>Float64,
    "oldmin"=>Float64,
    "old max"=>Float64,
    "oldmax"=>Float64,
    "line skip"=>Int,
    "lineskip"=>Int,
    "byte skip"=>Int,
    "byteskip"=>Int,
    # orientation
    "space dimension"=>Int,
    "space units"=>VTuple{QString},
    "space origin"=>PTuple{IntFloat},
    "space directions"=>VTuple{StringPTuple{Float64}},
    "measurement frame"=>VTuple{PTuple{Float64}},
    # per-axis
    "sizes"=>VTuple{Int},
    "spacings"=>VTuple{IntFloat},
    "thicknesses"=>VTuple{Float64},
    "axis mins"=>VTuple{IntFloat},
    "axismins"=>VTuple{IntFloat},
    "axis maxs"=>VTuple{IntFloat},
    "axismaxs"=>VTuple{IntFloat},
    "centers"=>VTuple{String},
    "centerings"=>VTuple{String},
    "labels"=>VTuple{QString},
    "units"=>VTuple{QString},
    "kinds"=>VTuple{String},
)

const fieldorder = ["content", "type", "dimension", "space", "space dimension",
                    "sizes", "spacings", "space directions", "kinds",
                    "centers", "centerings", "thickness",
                    "axis mins", "axismins", "axis maxs", "axismaxs",
                    "labels", "units",
                    "min", "max", "old min", "oldmin", "old max", "oldmax",
                    "block size", "blocksize", "endian", "encoding",
                    "space units", "space origin", "measurement frame",
                    "line skip", "lineskip", "byte skip", "byteskip",
                    "sample units", "sampleunits",
                    "data file", "datafile"]

const per_axis = [# orientation-related fields
                  "space directions",
                  # other per-axis fields
                  "sizes", "spacings", "thicknesses", "axis mins", "axismins",
                  "axis maxs", "axismaxs", "centers", "centerings", "labels",
                  "units", "kinds"]

const per_spacedim = ["space units", "space origin", "measurement frame"]

# version >= n is required if it has any fields in version_reqs[n]
const version_reqs = ([],
                      [],   # key/value tested separately,
                      ["kinds"],
                      ["thicknesses", "sample units", "space",
                       "space dimension", "space directions",
                       "space origin", "space units",
                       "data file", "datafile"],
                      ["measurement frame"],
)

function myendian()
    if ENDIAN_BOM == 0x04030201
        return "little"
    elseif ENDIAN_BOM == 0x01020304
        return "big"
    end
end

# Don't extend FileIO.load
# Set mode to "r+" if you want to be able to modify values in the
# image and have them update in the disk file
function load(f::File{format"NRRD"}, args...; mode="r", mmap=:auto)
    open(f, mode) do io
        skipmagic(io)
        load(io, args...; mode=mode, mmap=mmap)
    end
end

function load(io::Stream{format"NRRD"}, Tuser::Type=Any; mode="r", mmap=:auto)
    # Assemble all the information about the array we're about to
    # read: element type, size, and the "meaning" of axes
    version, header, keyvals, comments = parse_header(io)
    Traw, need_bswap = raw_eltype(header)
    szraw = (header["sizes"]...,)  # "sizes" may change in outer_eltype!, grab it now
    T, nd, perm = outer_eltype!(header, Traw, Tuser)
    axs = get_axes(header, nd)
    sz = get_size(axs)

    # Read the data
    iodata = find_datafile(io, header; mode=mode)
    compressed = in(header["encoding"], ("gzip", "gz"))
    if compressed
        iodata = Libz.ZlibInflateInputStream(iodata)
    end

    can_mmap = header["encoding"] == "raw"

    if mmap == true && (!can_mmap)
        error("Cannot use memory-mapped for reading a non-raw or bswapped file")
    end

    # Use memory-mapping for large files
    do_mmap = can_mmap && (prod(szraw) > 10^8) && (mmap == :auto)
    do_mmap |= can_mmap && (mmap == true)

    if !compressed
        szraw, sz = checked_size(Traw, szraw, sz, iodata)
    end

    if do_mmap
        # Recent Julia versions are picky about alignment
        pos = position(iodata)
        if pos % sizeof(Traw) != 0
            szcor = (szraw[1]*sizeof(Traw), szraw[2:end]...)
            A = reinterpret(Traw, Mmap.mmap(iodata, Array{UInt8,length(szraw)}, szcor, pos;
                          grow=false))
        else
            A = Mmap.mmap(iodata, Array{Traw,length(szraw)}, szraw, pos;
                          grow=false)
        end
        if need_bswap
            f = mode == "r+" ? (bswap, bswap) : bswap
            A = mappedarray(f, A)
        end
    elseif header["encoding"] == "raw" || in(header["encoding"], ("gzip", "gz"))
        A = read!(iodata, Array{Traw}(undef, szraw...))
        if need_bswap
            A = [bswap(a) for a in A]
        end
    else
        error("\"", header["encoding"], "\" encoding not supported.")
    end

    if perm == ()
        if T != eltype(A)
            A = need_bswap ? A = mappedarray(x->T(x), A) : reshape(reinterpret(T, A), sz)
        end
    else
        A = permuteddimsview(A, perm)
        if T<:Color
            A = colorview(T, A)
        end
    end

    isa(axs, Dims) ? A : AxisArray(A, axs)
end

function save(f::File{format"NRRD"}, img::AbstractArray; kwargs...)
    open(f, "w") do io
        write(io, magic(format"NRRD"))
        save(io, img; kwargs...)
    end
end

function save(io::Stream{format"NRRD"}, img::AbstractArray{T}; props::Dict = Dict{String,Any}(), keyvals=nothing, comments=nothing, kwargs...) where T
    axs = axisinfo(img)
    header = headerinfo(T, axs; kwargs...)
    header_eltype!(header, T)
    # copy fields from props to override those in header
    for (k, v) in props
        header[k] = v
    end
    v = version(header, !(keyvals==nothing || isempty(keyvals)))
    write_header(io, v, header, keyvals, comments)
    datafilename = get(props, "datafile", "")
    if isempty(datafilename)
        datafilename = get(props, "data file", "")
    end
    if isempty(datafilename)
        nrrd_write(io, img)
    else
        println(io.io, "data file: ", datafilename)
        if !get(props, "headeronly", false)
            open(datafilename, "w") do file
                nrrd_write(file, img)
            end
        end
    end
end

axisinfo(img) = axisinfo(HasAxes(img), img)
axisinfo(::HasAxes{true}, img) = AxisArrays.axes(img)
axisinfo(::HasAxes, img) = size(img)

### Interpreting header settings

"""
    arraytype!(header, version) -> T, axs, perm, need_bswap

Analyze the `header` dictionary to extract the element-type `T`, size
or axes information `axs`, the permutation `perm` (if any) that julia
should use for "wrapping" the read data, and a boolean `need_bswap`
indicating whether the data need to be byte-swapped (to account for
differences in endianness). `T` includes any color information (in
which case a dimension of the array will be "consumed"). `axs` will be
a Dims-tuple in simple cases, or an `Axes` tuple (from AxisArrays.jl)
if dimensions are labeled or have their spatial information
(pixelspacing, spacedirections, etc) specified. `perm` is the
permutation needed to move the color data to the first dimension, or
an empty tuple if no permutation is required.

This function may modify the `header` dictionary (the reason for the !
in the name), so make a copy first if necessary.
"""
function arraytype!(header, version)
    Traw, need_bswap = raw_eltype(header)
    T, nd, perm = outer_eltype!(header, Traw)
    axs = get_axes(header, nd)
    T, axs, perm, need_bswap
end

"""
    arraytype(filename)

Parse NRRD header and calls `arraytype!(header, version)`. See
`arraytype!` for information about the return values.
"""
function arraytype(filename)
    version, header, keyvals, comments = parse_header(filename)
    arraytype!(header, version)
end

function headerinfo(T, axs)
    header = Dict{String,Any}()
    Traw = raw_eltype(T)
    header["type"] = type2string(Traw)
    header["endian"] = myendian()
    header["encoding"] = "raw"
    if T <: Gray
        val = gray(oneunit(T))
        val = isa(val, FixedPoint) ? reinterpret(val) : val
        header["sample units"] = string("gray ", val)
    elseif T <: Union{RGB,RGBA}
        val = red(oneunit(T))
        val = isa(val, FixedPoint) ? reinterpret(val) : val
        valfmt = isa(val, Integer) ? 'd' : 'f'
        if T <: RGB
            colstr = "rgb"
            valfmtstr = "(%$valfmt,%$valfmt,%$valfmt)"
            vals = (val,val,val)
        else
            colstr = "rgba"
            valfmtstr = "(%$valfmt,%$valfmt,%$valfmt,%$valfmt)"
            vals = (val,val,val,val)
        end
        fmtstr = "%s $valfmtstr"
        x = (colstr, vals...)
        header["sample units"] = @eval @sprintf($fmtstr, $(x...))
    end
    # Do the axes information
    header["dimension"] = length(axs)
    if isa(axs, Base.Indices)
        axs = map(length, axs)
    end
    specifyorientation = false
    if isa(axs, Dims)
        header["sizes"] = [axs...]
    else
        # axs is an Axis-tuple
        header["sizes"] = [map(length, axs)...]
        axnames = map(ax->axisnames(ax)[1], axs)
        isspace = map(s->!startswith(string(s), "time"), axnames)
        if haskey(axes2space, axnames)
            header["space"] = axes2space[axnames]
            specifyorientation = true
        end
        header["kinds"] = [isspc ? "domain" : "time" for isspc in isspace]
        if !all(isdefaultname, axnames)
            header["labels"] = [string(s) for s in axnames]
        end
        rng = map(ax->axisvalues(ax)[1], axs)
        stepval = map(step, rng)
        unitstr = map(x->isa(x, Quantity) ? string(unit(x)) : "", stepval)
        spacing = map(x->isa(x, Quantity) ? ustrip(x) : x,        stepval)
        if !all(x->x=="", unitstr)
            header["units"] = [unitstr...]
        end
        if !all(x->x==1, spacing)
            if specifyorientation
                header["space directions"] = [ntuple(d->d==i ? spacing[i] : 0, length(spacing)) for i = 1:length(spacing)]
            else
                header["spacings"] = [Float64.(spacing)...]
            end
        end
        origin = map(x->isa(x, Quantity) ? ustrip(x) : x, map(first, rng))
        if specifyorientation && any(x->x!=0, origin)
            header["space origin"] = [origin[[isspace...]]...]
        end
    end
    # Adjust the axes for color
    if T <: Colorant && !(T <: AbstractGray)
        header["dimension"] = length(axs)+1
        pushfirst!(header["sizes"], length(T))
        if haskey(header, "spacings")
            pushfirst!(header["spacings"], NaN)
        end
        if haskey(header, "labels")
            pushfirst!(header["labels"], lowercase(string(T.name.name)))
        end
        if haskey(header, "units")
            pushfirst!(header["units"], "")
        end
        if !haskey(header, "kinds")
            header["kinds"] = ["domain" for d = 1:length(axs)]
        end
        if T <: RGB
            colkind = "RGB-color"
        elseif T <: HSV
            colkind = "HSV-color"
        elseif T <: XYZ
            colkind = "XYZ-color"
        elseif T <: RGBA
            colkind = "RGBA-color"
        else
            colkind = string(length(T), "-color")
        end
        pushfirst!(header["kinds"], colkind)
    end
    header
end

function version(header, has_keyvalue::Bool=false)
    for n = length(version_reqs):-1:1
        vr = version_reqs[n]
        for f in vr
            if haskey(header, f)
                return n
            end
        end
    end
    has_keyvalue ? 2 : 1
end

function isdefaultname(s::AbstractString)
    startswith(s, "dim_") || startswith(s, "space") ||
        startswith(s, "time") || startswith(s, "domain")
end
isdefaultname(s::Symbol) = isdefaultname(string(s))

"""
    raw_eltype(header) -> Traw, need_bswap
    raw_eltype(::Type{T}) -> Traw

Get the "basic" element type of the data, e.g., `UInt16` or
`Float32`.

This function does not try to determine whether the image is color
(`Traw` does not contain any color information), nor does it try to
interpret `Traw` as a `Normed` type.

See also: outer_eltype!, fixedtype.
"""
function raw_eltype(header)
    Traw = string2type[header["type"]]
    need_bswap = haskey(header, "endian") && header["endian"] != myendian() && sizeof(Traw) > 1
    Traw, need_bswap
end

raw_eltype(::Type{C}) where {C<:Colorant} = raw_eltype(eltype(C))
raw_eltype(::Type{T}) where {T<:FixedPoint}   = FixedPointNumbers.rawtype(T)
raw_eltype(::Type{T}) where {T} = T

"""
    fixedtype(Traw, header) -> Tu

Attempt to interpret type `Traw` in terms of FixedPoint numbers. The
interpretation depends on whether `header` has a "sample units" field
of the form:

    sample units: <colorspace> <whitepoint>

There must be a space between `colorspace` and `whitepoint`, and
`colorspace` must be one of "gray", "rgb", "rgba", "hsv", "xyz", or
their uppercase variants.  (Other than "gray", all of these are
supported "kinds" values. "gray" does not typically correspond to an
axis, which is why it isn't encoded in "kinds".) The presence of any
of these words indicates that the data represent an image rather than
some other kind of array.

Any other choices are ignored, as "sample units" can also be an
arbitrary string.

# Examples:

    # conventional uint8 grayscale
    sample units: gray 255

    # a 14-bit grayscale camera (numbers can be represented in hex format)
    sample units: gray 0x3fff

    # RGB encoded with float or double
    sample units: rgb (1.0,1.0,1.0)

    # RGB encoded with float or double, but using the scaling of uint8
    sample units: rgb (255.0,255.0,255.0)

    # conventional XYZ
    sample units: xyz (95.047,100.000,108.883)

    # HSV, hue measured in degrees
    sample units: hsv (360, 0, 1)

    # HSV, hue normalized to [0, 1]
    sample units: hsv (1, 0, 1)

If `Traw` cannot be interpreted as `Normed`, `Tu = Traw`.
"""
function fixedtype(::Type{Traw}, header) where Traw<:Unsigned
    # Note that "max" is not useful in this context.
    # See https://sourceforge.net/p/teem/bugs/14/
    if haskey(header, "sample units")
        su = header["sample units"]
    elseif haskey(header, "sampleunits")
        su = header["sampleunits"]
    else
        return Traw
    end
    idxsplit = findfirst(isequal(' '), su)
    idxsplit === nothing && return Traw
    cm, rest = lowercase(su[1:idxsplit-1]), strip(su[idxsplit+1:end])
    if cm ∈ ("gray", "rgb", "rgba", "xyz", "hsv")
        # It looks like a whitepoint definition
        if cm == "gray"
            val = numberparse(rest)
            return Gray{fixedtype_max(Traw, val)}
        elseif cm ∈ ("rgb", "rgba")
            s = nrrd_parse(PTuple{String}, rest)
            val = map(numberparse, s)
            if all(v->v==val[1], val)
                return fixedtype_max(Traw, val[1])
            end
        end
    end
    Traw
end
fixedtype(::Type{Traw}, header) where {Traw} = Traw

function fixedtype_max(::Type{Traw}, mx) where Traw<:Unsigned
    fmx = log2(mx+1)
    if round(fmx) == fmx
        return Normed{Traw,round(Int,fmx)}
    end
    Traw
end

header_eltype!(header, ::Type) = header
function header_eltype!(header, ::Type{T}) where T<:FixedPoint
    header["sample units"] = string("gray ", reinterpret(one(T)))
    header
end
function header_eltype!(header, ::Type{C}) where C<:Colorant
    _header_eltype!(header, C, eltype(C))
    header
end
function _header_eltype!(header, ::Type{C}, ::Type{T}) where {C<:AbstractGray,T<:FixedPoint}
    header["sample units"] = string("gray ", reinterpret(one(T)))
end
function _header_eltype!(header, ::Type{C}, ::Type{T}) where {C<:AbstractGray,T}
    header["sample units"] = string("gray ", one(T))
end
function _header_eltype!(header, ::Type{C}, ::Type{T}) where {C<:AbstractRGB,T<:FixedPoint}
    o = reinterpret(one(T))
    header["sample units"] = "rgb ($o,$o,$o)"
end
function _header_eltype!(header, ::Type{C}, ::Type{T}) where {C<:AbstractRGB,T}
    o = one(T)
    header["sample units"] = "rgb ($o,$o,$o)"
end
function _header_eltype!(header, ::Type{C}, ::Type) where C<:XYZ
    header["sample units"] = "xyz (95.047,100.000,108.883)"
end
function _header_eltype!(header, ::Type{C}, ::Type) where C<:HSV
    header["sample units"] = "hsv (360, 0, 1)"
end

"""
    colorant_eltype(C, T) -> Tc

Return a valid "inner" element type `Tc` for colorant type `C`. When
`T` != `Tc`, values must be "converted" before they can be interpreted
as type `C`.
"""
colorant_eltype(::Type{C}, ::Type{T}) where {C<:Colorant, T<:AbstractFloat} = C{T}
colorant_eltype(::Type{C}, ::Type{T}) where {C<:Colorant, T} = C{Float32}

"""
    UnknownColor{T,N}

An unknown Color. This type gets returned when one of the "kind"
settings is "3-color" or "4-color".
"""
struct UnknownColor{T,N} <: Color{T,N}
    col::NTuple{N,T}
end

"""
    T, nd, perm = outer_eltype!(header, Traw)

Extract the julia array `eltype` `T`, the number of dimensions `nd`
**excluding** color/complex/vector/matrix element data, and any
permutation needed to put the eltype dimension first. Any dimensions
in the header corresponding to color (or if "kind" is set to one of
the vector types) will be "consumed" upon exit. `Traw` is the
element type as determined by `raw_eltype`.

See also: raw_eltype.
"""
function outer_eltype!(header, Traw, Tuser=Any)
    nd = header["dimension"]
    sz = header["sizes"]
    length(sz) == nd || error("parsing of sizes: $(header["sizes"]) is inconsistent with $nd dimensions")
    perm = ()
    T = Traw
    if !(T <: Tuser)
        if Tuser <: Union{AbstractRGB,ColorTypes.TransparentRGB}
            if T <: Unsigned
                T = Normed{T,8*sizeof(T)}
            end
            T = ccolor(Tuser, base_colorant_type(Tuser){T})
            length(T) == sz[1] || error("first dimension of size $(sz[1]), expected $(length(T))")
            # We've just handled the first dimension, so delete the per-axis data for it
            nd -= 1
            for fn in per_axis
                if haskey(header, fn)
                    deleteat!(header[fn], 1)
                end
            end
        else
            T = Tuser
        end
    end
    if haskey(header, "kinds")
        kinds = header["kinds"]
        length(kinds) == nd || error("parsing of kinds: $(header["kinds"]) is inconsistent with $nd dimensions")
        for i = 1:nd
            k = kinds[i]
            if k == "RGB-color"
                chksize(sz[i], 3)
                Tu = fixedtype(Traw, header)
                Tu = Tu == UInt8 ? N0f8 : (Tu == UInt16 ? N0f16 : Tu)
                T = RGB{Tu}
            elseif k == "HSV-color"
                chksize(sz[i], 3)
                T = colorant_eltype(HSV, fixedtype(Traw, header))
            elseif k == "XYZ-color"
                chksize(sz[i], 3)
                T = colorant_eltype(XYZ, fixedtype(Traw, header))
            elseif k == "RGBA-color"
                chksize(sz[i], 4)
                Tu = fixedtype(Traw, header)
                T = RGBA{Tu}
            elseif k == "3-color"
                chksize(sz[i], 3)
                T = UnknownColor{Traw,3}
            elseif k == "4-color"
                chksize(sz[i], 4)
                T = UnknownColor{Traw,4}
            elseif k == "complex"
                chksize(sz[i], 2)
                T = Complex{Traw}
            elseif k == "quaternion"
                checksize(sz[i], 2)
                T = Quaternion{Traw}
            elseif k == "2-vector"
                chksize(sz[i], 2)
                T = SVector{2,Traw}
            elseif k == "3-vector" || k == "3-gradient" || k == "3-normal"
                chksize(sz[i], 3)
                T = SVector{3,Traw}
            elseif k == "4-vector"
                chksize(sz[i], 4)
                T = SVector{4,Traw}
            elseif k == "2D-matrix"
                chksize(sz[i], 4)
                T = SMatrix{2,2,Traw}
            elseif k == "3D-matrix"
                chksize(sz[i], 9)
                T = SMatrix{3,3,Traw}
            end
            if T != Traw
                # We've handled this dimension, adjust the dimensionality
                if i > 1
                    perm = (i, setdiff(1:nd, i)...)
                end
                nd -= 1
                for fn in per_axis
                    if haskey(header, fn)
                        deleteat!(header[fn], i)
                    end
                end
                break
            end
        end
    end
    if T == Traw && nd == header["dimension"]
        T = fixedtype(Traw, header)
    end
    T, nd, perm
end

function get_axes(header, nd)
    # Validate the per-axis fields...
    for f in per_axis
        if haskey(header, f)
            length(header[f]) == nd || error("expected $nd (remaining) dimensions in field $f, got $(header[f])")
        end
    end
    # ...and the orientation fields
    if haskey(header, "space")
        sd = space2axes[lowercase(header["space"])][1]
    else
        sd = get(header, "space dimension", nd)
    end
    for f in per_spacedim
        if haskey(header, f)
            length(header[f]) == sd || error("expected $sd dimensions in field $f, got $(header[f])")
        end
    end
    ## Test whether we have fields that require an AxisArray
    axes_fields = setdiff(per_axis, ["sizes", "kinds", "centers", "centerings"])
    need_axes = haskey(header, "space")
    for f in axes_fields
        need_axes |= haskey(header, f)
    end
    # also check for a time axis
    istime = falses(nd)
    isspace = falses(nd)
    if haskey(header, "kinds")
        kinds = header["kinds"]
        istime = map(x->x=="time", kinds)
        isspace = map(x->x=="space" || x=="domain", kinds)
        need_axes |= any(istime)
    elseif haskey(header, "space dimension")
        length(isspace) == sd || error("confused, have \"space dimension\" $sd but got $nd dimensions: $header")
        fill!(isspace, true)
    elseif haskey(header, "space")
        nd == sd || error("confused, have \"space\" but can't tell which axes are spatial: $header")
        fill!(isspace, true)
    end
    sz = header["sizes"]
    if !need_axes
        return (sz...,)
    end
    ## Axis names
    axnames = [string("dim_", d) for d = 1:nd]
    if haskey(header, "labels")
        axnames = header["labels"]
        if any(istime)
            for idx in findall(istime)
                labelt = axnames[idx]
                if !istimeaxis(Axis{Symbol(labelt)})
                    @warn("label $labelt is not defined as a time axis, define it with `@traitimpl TimeAxis{Axis{:$labelt}}` (see ImageAxes for more information)")
                end
            end
        end
    elseif haskey(header, "space")
        spcnames = map(string, space2axes[lowercase(header["space"])][2])
        copy_space!(axnames, spcnames, isspace)
    elseif haskey(header, "kinds")
        # Give axes default names based on their kind: space1, space2, etc.
        kindcounter = Dict{String,Int}()
        for (i,k) in enumerate(header["kinds"])
            if !haskey(kindcounter, k)
                kindcounter[k] = 0
            end
            n = kindcounter[k]+1
            kindcounter[k] = n
            # append numeric label if more than 1
            axnames[i] = n > 1 ? string(k, '_', n) : k
        end
        # ...but if there was more than 1, go back and append "1" at
        # the end of the first one
        for (k, n) in kindcounter
            if n > 1
                idx = findall(x->x==k, axnames)
                axnames[idx] .= string(k, '_', 1)
            end
        end
    end
    ## Ranges
    # Can a "space directions" field be encoded as a "spacings" field?
    if haskey(header, "space directions")
        can_convert, spc = spacings(header["space directions"])
        if all(isnan, spc)
            # they were all "none"
            delete!(header, "space directions")
            can_convert = false
        end
        if can_convert
            !haskey(header, "spacings") || error("cannot have both \"space directions\" and \"spacings\"")
            header["spacings"] = spc
            delete!(header, "space directions")
        end
    end
    if haskey(header, "axis mins") || haskey(header, "axismins")
        axmin = haskey(header, "axis mins") ? header["axis mins"] : header["axismins"]
        axmax = haskey(header, "axis maxs") ? header["axis maxs"] : header["axismaxs"]
        rng = Any[startstoplen(axmin[d], axmax[d], sz[d]) for d = 1:nd]
    else
        startval = Any[zeros(Int, nd)...]
        stepval = Any[ones(Int, nd)...]
        if haskey(header, "space origin")
            so = header["space origin"]
            copy_space!(startval, header["space origin"], isspace)
        end
        if haskey(header, "spacings")
            copy_not_nan!(stepval, header["spacings"])
        end
        rng = Any[startsteplen(startval[d], stepval[d], sz[d]) for d = 1:nd]
    end
    # These are more laborious than they should be because of issues
    # with rounding and ranges in julia-0.5.
    if haskey(header, "units")
        ua = [unit_string_dict[x] for x in header["units"]]
        rng = Any[startstepstop(first(r)*u, step(r)*u, last(r)*u) for (r,u) in zip(rng,ua)]
    elseif haskey(header, "space units")
        us = [get(unit_string_dict, x, 1) for x in header["space units"]]
        ua = fill!(Array{Any}(undef, nd), 1)
        copy_space!(ua, us, isspace)
        for d = 1:nd
            r, u = rng[d], ua[d]
            newr = startstepstop(first(r)*u, step(r)*u, last(r)*u)
            rng[d] = newr
        end
    end
    l = map(safelength, rng)
    l == sz || error("expect axis lengths to match array sizes, got $((rng...,)) of length $l and $sz")
    # Return the axes
    ntuple(d->Axis{Symbol(axnames[d])}(rng[d]), nd)
end

get_size(sz::Dims) = sz
get_size(axs) = map(safelength, axs)

### Parsing

"""
    parse_header(io) -> version, header, keyvals, comments

Parse the NRRD header (the top of the .nrrd or the separate .nhdr
file). `io` should be positioned just after the initial "NRRD" in the
file. This reads up to and including the first blank line of the file,
so at the end (if this is a properly-formatted NRRD file) `io` is
positioned at the first byte of the data (if present).

Outputs:
- `version` is a 4-character string, e.g., "0002", giving the NRRD version of the header.
- `header` is a `Dict{String,Any}` of `field=>setting` pairs (the settings are parsed, they are not necessarily strings)
- `keyvals` is a `Dict{String,String}` containing `key=>value` pairs (NRRD0002 or higher, lines like key:=value; many NRRD files do not contain any of these)
- `comments` is an array containing lines of the header that began with `#` (but with the `#` and leading whitespace stripped out)

See also: write_header.
"""
function parse_header(io)
    version = ascii(String(read!(io, Vector{UInt8}(undef, 4))))
    skipchars(isspace, io)
    header = Dict{String, Any}()
    keyvals = Dict{String, String}()
    comments = String[]
    # Read until we encounter a blank line, which is the separator
    # between the header and data
    line = strip(readline(io))
    while !isempty(line)
        if line[1] != '#'
            idx = findfirst(isequal(':'), line)
            idx === nothing && error("no colon found in $line")
            key, value = line[1:idx-1], line[idx+1:end]
            if !isempty(value) && value[1] == '='
                # This is a NRRD key/value pair, insert into keyvals
                keyvals[key] = value[2:end]
            else
                lkey = lowercase(key)
                T = get(parse_type, lkey, String)
                header[lkey] = nrrd_parse(T, strip(value))
            end
        else
            cmt = strip(lstrip(line, ['#', ' ']))
            if !isempty(cmt)
                push!(comments, cmt)
            end
        end
        line = strip(readline(io))
    end
    version, header, keyvals, comments
end

parse_header(s::Stream{format"NRRD"}) = parse_header(stream(s))

function parse_header(filename::AbstractString)
    f = File{format"NRRD"}(filename)
    open(f) do io
        skipmagic(io)
        parse_header(io)
    end
end

"""
    write_header(io, version, header, [keyvals, [comments]])

Write an NRRD header, as the top of the .nrrd or the separate .nhdr
file. `io` should be positioned just after the initial "NRRD" in the
file. This writes the header and a blank line, so that at the end `io`
is positioned at the first byte of the data (if present).

Note that if you're writing a header for a "detached" data file
(separate .nhdr and .raw files), `header` should contain a "data file"
or "datafile" field storing the name of the .raw file.

Inputs:
- `version` is a 4-character string, e.g., "0002", giving the NRRD version of the header, or a integer corresponding to a recognized NRRD header version number
- `header` is a `Dict{String,Any}` of `field=>setting` pairs (as returned by `parse_header`)
- `keyvals` is a `Dict{String,String}` containing `key=>value` pairs (NRRD0002 or higher, lines like key:=value; many NRRD files do not contain any of these)
- `comments` is an array containing lines of the header that began with `#` (but with the `#` and leading whitespace stripped out)

See also: parse_header.
"""
function write_header(io::IO, version, header, keyvals=nothing, comments=nothing)
    writeversionstr(io, version)

    for fn in fieldorder
        if haskey(header, fn)
            print(io, fn, ": ")
            T = get(parse_type, fn, String)
            nrrd_format(io, T, header[fn])
            println(io, "")
        end
    end
    if keyvals != nothing
        for (k,v) in keyvals
            println(io, k, ":=", v)
        end
    end
    if comments != nothing
        for c in comments
            println(io, "# ", c)
        end
    end
    println(io)
    nothing
end
write_header(s::Stream{format"NRRD"}, args...) = write_header(stream(s), args...)

nrrd_parse(::Type{T}, str) where {T} = parse(T, str)  # fallback

function nrrd_parse(::Type{IntFloat}, str)
    x = parse(Float64, str)
    x == round(x) ? Int(x) : x
end

nrrd_parse(::Type{String}, str) = str

function nrrd_parse(::Type{QString}, str)
    str[1] == '"' && str[end] == '"' || error("$str does not start and end with double-quote")
    str[2:end-1]
end

function nrrd_parse(::Type{VTuple{T}}, s::AbstractString) where T
    ss = split(s)  # r"[ ,;]")
    v = Vector{alloctype(T)}(undef, length(ss))
    for i = 1:length(ss)
        v[i] = nrrd_parse(T, ss[i])
    end
    return v
end

function nrrd_parse(::Type{VTuple{QString}}, s::AbstractString)
    str = nrrd_parse(QString, s)  # to strip the first and last ""
    split(str, r"\" +\"")
end

function nrrd_parse(::Type{PTuple{T}}, s::AbstractString) where T
    s[1] == '(' && s[end] == ')' || error("$s should begin and end parentheses")
    str = s[2:end-1]
    ss = split(str, ',', keepempty=false)
    v = Vector{alloctype(T)}(undef, length(ss))
    for i = 1:length(ss)
        v[i] = nrrd_parse(T, ss[i])
    end
    return v
end

function nrrd_parse(::Type{StringPTuple{T}}, s::AbstractString) where T
    if s[1] == '(' && s[end] == ')'
        return nrrd_parse(PTuple{T}, s)
    end
    s
end

writeversionstr(io, v::Integer) = @printf(io, "%04d\n", v)
function writeversionstr(io, v::AbstractString)
    length(v) == 4 || error("version string must have 4 characters, got $v")
    println(io, v)
end

nrrd_format(io, ::Type{T}, x) where {T} = print(io, x)

function nrrd_format(io, ::Type{QString}, str)
    print(io, '"', str, '"')
end

function nrrd_format(io, ::Type{VTuple{T}}, container) where T
    isfirst = true
    for x in container
        if isfirst
            isfirst = false
        else
            print(io, ' ')
        end
        nrrd_format(io, T, x)
    end
    nothing
end

function nrrd_format(io, ::Type{PTuple{T}}, container) where T
    print(io, '(')
    isfirst = true
    for x in container
        if isfirst
            isfirst = false
        else
            print(io, ',')
        end
        nrrd_format(io, T, x)
    end
    print(io, ')')
end

nrrd_format(io, ::Type{StringPTuple{T}}, s::AbstractString) where {T} = nrrd_format(io, String, s)
nrrd_format(io, ::Type{StringPTuple{T}}, container) where {T} = nrrd_format(io, PTuple{T}, container)

alloctype(::Type{T}) where {T} = T
alloctype(::Type{PTuple{T}}) where {T} = Any
alloctype(::Type{StringPTuple{T}}) where {T} = Any # Union{String,Vector{T}}
alloctype(::Type{IntFloat}) = Union{Int,Float64}

### Utilities

function chksize(sz, sztarget)
    sz == sztarget || error("dimension size should be $sztarget, got $sz")
    nothing
end

function startstoplen(min::Integer, max::Integer, len::Integer)
    if max-min == len-1
        return Int(min):Int(max)
    end
    stepval = (max-min)/(len-1)
    if stepval ≈ round(stepval)
        return Int(min):round(Int,stepval):Int(max)
    end
    startstopsteplen(min, max, stepval, len)
end

function startstoplen(min, max, len)
    stepval = (max-min)/(len-1)
    startstopsteplen(min, max, stepval, len)
end

function startsteplen(min, stepval, len)
    max = min+stepval*(len-1)
    startstoplen(min, max, len)
end

function startstepstop(min, stepval, max)
    n = (max-min)/stepval + 1
    len = floor(Int, n+10*eps(n))
    startstopsteplen(min, max, stepval, len)
end

function startstopsteplen(min::T, max::T, stepval::T, len) where T<:Integer
    stepval == 1 ? (min:max) : (min:stepval:max)
end
function startstopsteplen(min::T, max::T, stepval::T, len) where T<:AbstractFloat
    imin, imax, istepval = round(Int, min), round(Int, max), round(Int, stepval)
    if min ≈ imin && max ≈ imax && stepval ≈ istepval
        return startstopsteplen(imin, imax, istepval, len)
    end
    r = range(min, stop=max, length=len)
    abs(step(r)/stepval - 1) < 0.01 || error("wanted length $len from $min to $max, got $r of length $(length(r))")
    r
end
function startstopsteplen(min::T, max::T, stepval::T, len) where T
    r = range(min, stop=max, length=len)
    abs(step(r)/stepval - 1) < 0.01 || error("wanted length $len from $min to $max, got $r of length $(length(r))")
    r
end
startstopsteplen(min, max, stepval, len) = startstopsteplen(promote(min, max, stepval)..., len)

function copy_space!(dest, src, isspace)
    length(src) == sum(isspace) || throw(DimensionMismatch("src length $(length(src)) disagrees with isspace $isspace"))
    length(dest) == length(isspace) || throw(DimensionMismatch("dest length $(length(dest)) does not agree with isspace of length $(length(isspace))"))
    d = 1
    for ds = 1:length(src)
        while !isspace[d]
            d += 1
        end
        dest[d] = src[ds]
        d += 1
    end
    dest
end

function copy_not_nan!(dest, src)
    nd = length(dest)
    length(src) == nd || throw(DimensionMismatch("dest length $nd does not match source length $(length(src))"))
    for d = 1:nd
        s = src[d]
        if !isa(s, AbstractFloat) || !isnan(s)
            dest[d] = s
        end
    end
    dest
end

function spacings(spacedirections)
    spacings = zeros(length(spacedirections))
    can_convert = true
    d = 1
    for (i, v) in enumerate(spacedirections)
        if v == "none"
            spacings[i] = NaN
        else
            n = 0
            for (j,x) in enumerate(v)
                if x != 0
                    if j == d   # make certain they are in order
                        spacings[i] = x
                        d += 1
                    else
                        can_convert = false
                    end
                    n += 1
                end
            end
            can_convert &= n < 2
        end
    end
    can_convert, spacings
end

### File-related functions

# Adjust the array size if the file is not big enough for reading to succeed
function checked_size(Traw, szraw, sz, iodata)
    cpos = position(iodata)
    datalen = div(filesize(stat(iodata)) - cpos, sizeof(Traw))
    if datalen < prod(szraw)
        if szraw != sz
            # If we've dropped a dimension due to "absorbing"
            # color, etc, don't try to figure out how to correct
            # it, just punt to the user
            error("data length $datalen too small for size $sz array of $T")
        end
        # If the data are smaller than the header suggests, read in as
        # much "complete" data as are available
        sznew = [szraw...]
        strds = [1;cumprod(sznew)]
        k = length(sznew)
        sznew[k] = div(datalen, strds[k])
        while sznew[k] == 0 && k > 1
            pop!(sznew)
            k -= 1
            sznew[k] = div(datalen, strds[k])
        end
        tsznew = (sznew...,)
        @warn("header indicates an array size $szraw, but the file size is consistent with at most $tsznew")
        szraw = sz = tsznew
    end
    szraw, sz
end

function stream2name(s::IO)
    name = s.name
    if !startswith(name, "<file ")
        error("io name ", name, " doesn't fit expected pattern")
    end
    name[7:end-1]
end

function find_datafile(iodata::IOStream, header; mode="r", path=nothing)
    if haskey(header, "data file") || haskey(header, "datafile")
        # Separate header and data files
        fdata = haskey(header, "data file") ? header["data file"] : header["datafile"]
        if path == nothing
            # Check current directory first
            if isfile(fdata)
                iodata = open(fdata)
            else
                path = dirname(stream2name(iodata))
                iodata = open(joinpath(path, fdata), mode)
            end
        else
            iodata = open(joinpath(path, fdata), mode)
        end
    end
    iodata
end

function find_datafile(s::Stream{format"NRRD"}, header; mode="r")
    find_datafile(stream(s), header; mode=mode)
end

nrrd_write(io, A::AxisArray) = nrrd_write(io, A.data)
nrrd_write(io, A::AbstractArray) = nrrd_write_elty(io, A)
nrrd_write_elty(io, A::AbstractArray{C}) where {C<:Colorant} = nrrd_write_elty(io, channelview(A))
nrrd_write_elty(io, A::AbstractArray{T}) where {T<:FixedPoint} = nrrd_write_elty(io, rawview(A))
nrrd_write_elty(io, A::AbstractArray) = write(io, A)

function numberparse(str)
    lstr = lowercase(str)
    startswith(lstr, "0x") && return parse(UInt, lstr)   # hex
    val = parse(Float64, lstr)
end

# Avoids roundoff error in Base's implementation of length
function safelength(r::StepRange)
    n = ((last(r)-first(r))/step(r)) + 1
    floor(Int, n+10*eps(n))
end
safelength(r) = length(r)

end
