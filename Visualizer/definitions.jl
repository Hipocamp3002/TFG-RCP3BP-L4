using StaticArrays

struct intersection
    u::Union{Nothing,SVector{4,Float64}}
    t::Float64
end

muDirectory(mu_val) = "./data/"*currSection*"/"*string(mu_val)*".jld2"
