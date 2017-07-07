# NRRD
[![NRRD](http://pkg.julialang.org/badges/NRRD_0.4.svg)](http://pkg.julialang.org/?pkg=NRRD)
[![NRRD](http://pkg.julialang.org/badges/NRRD_0.5.svg)](http://pkg.julialang.org/?pkg=NRRD)
[![Build Status](https://travis-ci.org/JuliaIO/NRRD.jl.svg?branch=master)](https://travis-ci.org/JuliaIO/NRRD.jl)
[![codecov.io](http://codecov.io/github/JuliaIO/NRRD.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaIO/NRRD.jl?branch=master)

Package for reading NRRD files.
Implements the  [FileIO](https://github.com/JuliaIO/FileIO.jl) interface.

Originally located in [Images.jl](https://github.com/timholy/Images.jl)

## Writing plain NRRD headers

Normal usage is as easy as

```julia
img = load("myfile.nrrd")
img = load("myfile.nhdr")
save("myotherfile.nrrd", img)
```

However, if you already have a raw binary representing the "data
file", the FileIO interface isn't sufficently flexible for writing
just the header. Assuming you want to save "rich" axis information, a
low-level approach using AxisArrays is the following:

```julia
using NRRD, FileIO, FixedPointNumbers, AxisArrays, Unitful
using Unitful: μm, s

# For a 480x640x200 image with time as the third axis,
# assuming a pixel spacing of 0.25μm and a framerate of 8fps
axy = Axis{:y}((1:480)*0.25μm)
axx = Axis{:x}((1:640)*0.25μm)
axt = Axis{:time}((1:200)*0.125s)

header = NRRD.headerinfo(N0f16, (axy, axx, axt))  # assuming N0f16 data
header["datafile"] = "mydata.raw"

open("mydata.nhdr", "w") do io
    write(io, magic(format"NRRD"))
    NRRD.write_header(io, "0004", header)
end
```
