using FileIO, FixedPointNumbers, Colors, ImageCore, Base.Test

workdir = joinpath(tempdir(), "ImagesNRRD")
isdir(workdir) || mkdir(workdir)
const teembase = "http://teem.sourceforge.net/nrrd/files"
const base3D   = "http://www.tc18.org/code_data_set/3D_greyscale"

function getfile(name, base=teembase)
    file = joinpath(workdir, name)
    if !isfile(file)
        file = download(joinpath(base, name), file)
    end
    file
end

@testset "readremote" begin
    fool_file = getfile("fool.nrrd")
    imgu = load(fool_file)
    @test eltype(imgu) == UInt8
    @test size(imgu) == (128,128)
    @test imgu[120,5] < imgu[5,5] < imgu[120,120]  # a simple orientation test
    foolf_file = getfile("foolf.nrrd")
    imgf = load(foolf_file)
    @test eltype(imgf) == Float32
    @test size(imgf) == (128,128)
    @test imgf[120,5] < imgf[5,5] < imgf[120,120]
    @test maximum(abs(imgf-imgu))/maximum(imgu) <= 1/200
    foolc_file = getfile("foolc.nrrd")
    imgc = load(foolc_file, RGB)
    @test eltype(imgc) == RGB{N0f8}
    @test size(imgc) == (128,128)
    imgg = convert(Array{Gray}, imgc)
    @test mean(abs(Int.(rawview(channelview(imgg))) - Int.(imgu))) < 12
    aneurism_rawfile = getfile("aneurism.raw.gz", base3D)
    aneurism_file = joinpath(workdir, "aneurism.nhdr")
    aneurism_rawfile_path = joinpath(workdir, "aneurism.raw.gz")
    open(aneurism_file, "w") do io
        write(io, """
NRRD0001
content: aneurism
# Courtesy of Philips Research, Hamburg, Germany
dimension: 3
type: uchar
sizes: 256 256 256
spacings: 1 1 1
data file: $aneurism_rawfile_path
encoding: gzip
""")
    end
    imga = load(aneurism_file)
    @test eltype(imga) == UInt8
    @test size(imga) == (256,256,256)
    @test sdims(imga) == 3
    @test pixelspacing(imga) == (1.0,1.0,1.0)
end

gc()  # to close any mmapped files
rm(workdir, recursive=true)

nothing
