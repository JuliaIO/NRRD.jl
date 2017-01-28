using FileIO, NRRD, FixedPointNumbers, Colors

Alast = nothing
for (A, name) in ((rand(UInt8, 3, 5), "uint8"),
                  (reinterpret(N0f8, rand(UInt8, 3, 5)), "u8"),
                  (rand(Gray{N0f8}, 3, 5), "gray_u8"),
                  (rand(RGB{N0f8}, 3, 5), "rgb_u8"),
                  (rand(Gray{N0f16}, 3, 5), "gray_u16"),
                  (reinterpret(Gray{N4f12}, rand(0x0000:0x0fff, 3, 5)), "gray_u12"))
    fn = joinpath(dirname(@__FILE__), "io", "eltype_$name.nrrd")
    save(fn, A)
    Alast = A
end

# Now manually tweak for the "hex" version
fn = joinpath(dirname(@__FILE__), "io", "eltype_gray_u12_hex.nrrd")
save(fn, Alast, props=Dict("sample units" => "gray 0x0fff"))
