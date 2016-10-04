module NRRD

# Packages needed to return the possible range of element types
using FixedPointNumbers, Colors, StaticArrays, Quaternions
# Other packages
using AxisArrays, ImageAxes, Unitful, MappedArrays, Ranges
using FileIO
import Libz
import FixedPointNumbers

typedict = Dict(
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

spacedimdict = Dict(
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

# We put these in a dict so that we don't eval untrusted
# strings. Please submit PRs to add to this list if you need
# additional unit support.
unit_string_dict = Dict("m" => u"m", "mm" => u"mm", "s" => u"s",
                        "um" => u"μm", "μm" => u"μm")

immutable QString end                 # string with quotes around it: "mm"
typealias VTuple{T} Tuple{Vararg{T}}  # space-delimited tuple: 80 150
immutable PTuple{T} end               # parenthesis-delimited tuple: (80,150)
immutable StringPTuple{T} end         # string or PTuple{T}

typealias IntFloat Union{Int,Float64}

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
    "sizes"=>Dims,
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

const per_axis = [# orientation-related fields
                  "space directions",
                  # other per-axis fields
                  "sizes", "spacings", "thicknesses", "axis mins", "axismins",
                  "axis maxs", "axismaxs", "centers", "centerings", "labels",
                  "units", "kinds"]

const per_spacedim = ["space units", "space origin", "measurement frame"]

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
function load(f::File{format"NRRD"}; mode="r", mmap=:auto)
    open(f, mode) do io
        skipmagic(io)
        load(io, mmap=mmap; mode=mode)
    end
end

function load(io::Stream{format"NRRD"}; mode="r", mmap=:auto)
    # Assemble all the information about the array we're about to
    # read: element type, size, and the "meaning" of axes
    version, header, keyvals, comments = parse_header(io)
    Traw, need_bswap = raw_eltype(header)
    szraw = (header["sizes"]...,)  # "sizes" may change in outer_eltype!, grab it now
    T, nd, perm = outer_eltype!(header, Traw)
    axs = get_axes(header, nd)
    sz = get_size(axs)

    # Read the data
    iodata = find_datafile(io, header; mode=mode)
    if in(header["encoding"], ("gzip", "gz"))
        iodata = Libz.ZlibInflateInputStream(iodata)
    end

    can_mmap = header["encoding"] == "raw"

    if mmap == true && (!can_mmap)
        error("Cannot use memory-mapped for reading a non-raw or bswapped file")
    end

    # Use memory-mapping for large files
    do_mmap = can_mmap && (prod(szraw) > 10^8) && (mmap == :auto)
    do_mmap |= can_mmap && (mmap == true)

    if do_mmap
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
            warn("header indicates an array size $szraw, but the file size is consistent with at most $tsznew")
            szraw = tsznew
        end
        A = Mmap.mmap(iodata, Array{Traw,length(szraw)}, szraw, cpos; grow=false)
        if need_bswap
            f = mode == "r+" ? (bswap, bswap) : bswap
            A = mappedarray(f, A)
        end
    elseif header["encoding"] == "raw" || in(header["encoding"], ("gzip", "gz"))
        A = read(iodata, T, szraw...)
        if need_bswap
            A = [bswap(a) for a in A]
        end
    else
        error("\"", header["encoding"], "\" encoding not supported.")
    end

    A = reinterpret(T, A, sz)
    isa(axs, Dims) ? A : AxisArray(A, axs)
end

function FileIO.save(f::File{format"NRRD"}, img::AbstractArray; props::Dict = Dict{String,Any}())
    sheader = stream(open(f, "w"))
    println(sheader, "NRRD0001")
    # Write the datatype
    T = get(props, "type", eltype(data(img)))
    if T<:AbstractFloat
        println(sheader, "type: ", (T == Float32) ? "float" :
                                   (T == Float64) ? "double" :
                                   (T == Float16) ? "float16" :
                                   error("Can't write type $T"))
    elseif eltype(T) <: FixedPointNumbers.UFixed
        # we don't want to write something like "type: ufixedbase{uint8,8}", so fix it
        T = FixedPointNumbers.rawtype(eltype(eltype(img)))
        println(sheader, "type: ", lowercase(string(T)))
    else
        println(sheader, "type: ", lowercase(string(T)))
    end
    # Extract size and kinds
    sz = get(props, "sizes", size(img))
    kinds = ["space" for i = 1:length(sz)]
    td = timedim(img)
    if td != 0
        kinds[td] = "time"
    end
    cd = colordim(img)
    if cd != 0
        if colorspace(img) == "RGB"
            kinds[cd] = "RGB-color"
        elseif size(img, cd) == 3
            kinds[cd] = "3-color"
        else
            kinds[cd] = "list"
        end
    end
    kinds = get(props, "kinds", kinds)
    # Write size and kinds
    println(sheader, "dimension: ", length(sz))
    print(sheader, "sizes:")
    for z in sz
        print(sheader, " ", z)
    end
    print(sheader, "\nkinds:")
    for k in kinds
        print(sheader, " ", k)
    end
    print(sheader, "\n")
    println(sheader, "encoding: ", get(props, "encoding", "raw"))
    println(sheader, "endian: ", get(props, "endian", ENDIAN_BOM == 0x04030201 ? "little" : "big"))
    ps = get(props, "pixelspacing", pixelspacing(img))
    ps = convert(Vector{Any}, ps)
    index = sort([cd, td])
    index[1] > 0 && insert!(ps, index[1], "nan")
    index[2] > 0 && insert!(ps, index[2], "nan")
    print(sheader, "spacings:")
    printunits = false
    for x in ps
        if isa(x, SIUnits.SIQuantity)
            printunits = true
            print(sheader, " ", x.val)
        else
            print(sheader, " ", x)
        end
    end
    print(sheader,"\n")
    if printunits
        print(sheader, "space units:")
        for x in ps
            if isa(x, SIUnits.SIQuantity)
                print(sheader," \"", strip(string(SIUnits.unit(x))), "\"")
            end
        end
        print(sheader, "\n")
    end
    datafilename = get(props, "datafile", "")
    if isempty(datafilename)
        datafilename = get(props, "data file", "")
    end
    if isempty(datafilename)
        write(sheader, "\n")
        write(sheader, data(img))
    else
        println(sheader, "data file: ", datafilename)
        if !get(props, "headeronly", false)
            open(datafilename, "w") do file
                write(file, data(img))
            end
        end
    end
    close(sheader)
end

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

"""
    raw_eltye(header) -> T, need_bswap

Get the "basic" element type of the data, e.g., `UInt16` or
`Float32`.

This function does not try to determine whether the image is color
(`T` does not contain any color information), nor does it try to
interpret `T` as a `UFixed` type.

See also: outer_eltype!, fixedtype.
"""
function raw_eltype(header)
    Traw = typedict[header["type"]]
    need_bswap = haskey(header, "endian") && header["endian"] != myendian() && sizeof(Traw) > 1
    Traw, need_bswap
end

# """
#     fixedtype(Traw, header) -> Tu

# Attempt to interpret type `Traw` as a `UFixed` type. The
# interpretation depends in part on the `"max"` field of `header`.

# If `Traw` cannot be interpreted as `UFixed`, `Tu = Traw`.
# """
# function fixedtype{Traw<:Unsigned}(::Type{Traw}, header)
#     f = 8*sizeof(Traw)
#     if haskey(header, "max")
#         mx = parse(Traw, header["max"])
#         fmx = log2(mx+1)
#         if fmx == round(fmx)
#             f = round(Int, fmx)
#         end
#     end
#     UFixed{Traw,f}
# end

# "max" is not useful in this context. See https://sourceforge.net/p/teem/bugs/14/
fixedtype(::Type{UInt8}, header) = UFixed8
fixedtype{Traw}(::Type{Traw}, header) = Traw

"""
    colorant_eltype(C, T) -> Tc

Return a valid "inner" element type `Tc` for colorant type `C`. When
`T` != `Tc`, values must be "converted" before they can be interpreted
as type `C`.
"""
colorant_eltype{C<:Colorant, T<:AbstractFloat}(::Type{C}, ::Type{T}) = C{T}
colorant_eltype{C<:Colorant, T}(::Type{C}, ::Type{T}) = C{Float32}

"""
    UnknownColor{T,N}

An unknown Color. This type gets returned when one of the "kind"
settings is "3-color" or "4-color".
"""
immutable UnknownColor{T,N} <: Color{T,N}
    col::NTuple{N,T}
end

"""
    outer_eltype!(header, Traw) -> T, nd, perm

Extract the julia array `eltype` `T`, the number of dimensions `nd`
**excluding** color/complex/vector/matrix element data, and any
permutation needed to put the eltype dimension first. Any dimensions
in the header corresponding to color (or if "kind" is set to one of
the vector types) will be "consumed" upon exit. `Traw` is the
element type as determined by `raw_eltype`.

See also: raw_eltype.
"""
function outer_eltype!(header, Traw)
    nd = header["dimension"]
    sz = header["sizes"]
    length(sz) == nd || error("parsing of sizes: $(header["sizes"]) is inconsistent with $nd dimensions")
    perm = ()
    T = Traw
    if haskey(header, "kinds")
        kinds = header["kinds"]
        length(kinds) == nd || error("parsing of kinds: $(header["kinds"]) is inconsistent with $nd dimensions")
        for i = 1:nd
            k = kinds[i]
            if k == "RGB-color"
                chksize(sz[i], 3)
                Tu = fixedtype(Traw, header)
                T = RGB{Tu}
            elseif k == "HSV-color"
                chksize(sz[i], 3)
                T = colorant_eltype(HSV, Traw)
            elseif k == "XYZ-color"
                chksize(sz[i], 3)
                T = colorant_eltype(XYZ, Traw)
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
        sd = spacedimdict[lowercase(header["space"])][1]
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
    need_axes = false
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
            for idx in find(istime)
                labelt = axnames[idx]
                if !istimeaxis(Axis{Symbol(labelt)})
                    warn("label $labelt is not defined as a time axis, define it with `@traitimpl TimeAxis{Axis{:$labelt}}` (see ImageAxes for more information)")
                end
            end
        end
    elseif haskey(header, "space")
        spcnames = map(string, spacedimdict[lowercase(header["space"])][2])
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
                idx = find(x->x==k, axnames)
                axnames[idx] = string(k, '_', 1)
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
        us = [unit_string_dict[x] for x in header["space units"]]
        ua = fill!(Array{Any}(nd), 1)
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
- `header` is a `Dict{String,String}` of `field=>setting` pairs (the settings are not parsed, they are pure strings)
- `keyvals` is a `Dict{String,String}` containing `key=>value` pairs (NRRD0002 or higher lines like key:=value; most NRRD files do not seem to contain any of these)
- `comments` is an array containing lines of the header that began with `#` (but with the `#` and leading whitespace stripped out)
"""
function parse_header(io)
    version = ascii(String(read(io, UInt8, 4)))
    skipchars(io,isspace)
    header = Dict{String, Any}()
    keyvals = Dict{String, String}()
    comments = String[]
    # Read until we encounter a blank line, which is the separator
    # between the header and data
    line = strip(readline(io))
    while !isempty(line)
        if line[1] != '#'
            key, value = split(line, ":")
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

nrrd_parse{T}(::Type{T}, str) = parse(T, str)  # fallback

function nrrd_parse(::Type{IntFloat}, str)
    x = parse(Float64, str)
    x == round(x) ? Int(x) : x
end

nrrd_parse(::Type{String}, str) = str

function nrrd_parse(::Type{QString}, str)
    str[1] == '"' && str[end] == '"' || error("$str does not start and end with double-quote")
    str[2:end-1]
end

function nrrd_parse{T}(::Type{VTuple{T}}, s::AbstractString)
    ss = split(s)  # r"[ ,;]")
    v = Array(alloctype(T), length(ss))
    for i = 1:length(ss)
        v[i] = nrrd_parse(T, ss[i])
    end
    return v
end

function nrrd_parse(::Type{VTuple{QString}}, s::AbstractString)
    str = nrrd_parse(QString, s)  # to strip the first and last "
    split(str, r"\" +\"")
end

function nrrd_parse{T}(::Type{PTuple{T}}, s::AbstractString)
    s[1] == '(' && s[end] == ')' || error("$s should begin and end parentheses")
    str = s[2:end-1]
    ss = split(str, ',', keep=false)
    v = Array(T, length(ss))
    for i = 1:length(ss)
        v[i] = nrrd_parse(T, ss[i])
    end
    return v
end

function nrrd_parse{T}(::Type{StringPTuple{T}}, s::AbstractString)
    if s[1] == '(' && s[end] == ')'
        return nrrd_parse(PTuple{T}, s)
    end
    s
end

alloctype{T}(::Type{T}) = T
alloctype{T}(::Type{StringPTuple{T}}) = Any # Union{String,Vector{T}}

function chksize(sz, sztarget)
    sz == sztarget || error("dimension size should be $sztarget, got $sz")
    nothing
end

### Utilities

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

function startstopsteplen{T<:Integer}(min::T, max::T, stepval::T, len)
    stepval == 1 ? (min:max) : (min:stepval:max)
end
function startstopsteplen{T}(min::T, max::T, stepval::T, len)
    r = Ranges.linspace(min, max, len)
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

# function parse_quantity(s::AbstractString, strict::Bool = true)
#     # Find the last character of the numeric component
#     m = match(r"[0-9\.\+-](?![0-9\.\+-])", s)
#     if m == nothing
#         error("AbstractString does not have a 'value unit' structure")
#     end
#     val = float64(s[1:m.offset])
#     ustr = strip(s[m.offset+1:end])
#     if isempty(ustr)
#         if strict
#             error("AbstractString does not have a 'value unit' structure")
#         else
#             return val
#         end
#     end
#     val * _unit_string_dict[ustr]
# end

# Avoids roundoff error in Base's implementation of length
function safelength(r::StepRange)
    n = ((last(r)-first(r))/step(r)) + 1
    floor(Int, n+10*eps(n))
end
safelength(r) = length(r)

end
