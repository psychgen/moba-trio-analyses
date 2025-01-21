
#author: Espen Eilertsen
#calculates GRM for trio-GCTA

using SnpArrays, CSV, LinearAlgebra, DataFrames, HDF5

function write_grm(G::AbstractMatrix, filename::AbstractString)
    g = ugrm2vec(G)
    nme = "grm_u"
    h5open(filename * ".h5", "w") do file
        write(file, nme, g)
    end
end

function ugrm2vec(M::AbstractMatrix{T}) where {T<:AbstractFloat}
    k = size(M, 1)
    p = Int(k * (k + 1) / 2)
    v = zeros(T, p)
    l = 1
    @inbounds for i in 1:k # cols
        @inbounds for j in 1:i # rows
            v[l] = M[j, i]
            l += 1
        end
    end
    v
end

function block_indices(size::Int, nblocks::Int)
    if nblocks > size
        throw(ArgumentError("Number of blocks can not be greater than the number of indices"))
    end
    T = promote_type(eltype(size), eltype(nblocks))
    blocksize = div(size, nblocks)
    rest = mod(size, nblocks)
    indices = Vector{UnitRange{T}}(undef, nblocks)
    # if size is not a multiple of nblocks, the first block gets the rest
    indices[1] = 1:blocksize + rest
    for i ∈ 2:nblocks
        first = last(indices[i - 1]) + 1
        indices[i] = first:first + blocksize - 1
    end
    indices
end

function grm_blocks(s::SnpArray, nblocks::Int)
    T = Float64
    npeople, nsnps = size(s)
    indices = block_indices(nsnps, nblocks)
    Φ = zeros(T, npeople, npeople)
    α = inv(2nsnps)

    # First block
    G = zeros(T, npeople, length(indices[1]))
    @views Base.copyto!(G, s[:, indices[1]], model=ADDITIVE_MODEL, impute=true, center=true, scale=true)
    BLAS.syrk!('U', 'N', α, G, one(T), Φ)
    println("Block: 1, indices: $(indices[1])")

    # Rest
    if(length(indices[2]) < length(indices[1]))
        G = zeros(T, npeople, length(indices[2]))
    end
    for i ∈ 2:nblocks
        @views Base.copyto!(G, s[:, indices[i]], model=ADDITIVE_MODEL, impute=true, center=true, scale=true)
        BLAS.syrk!('U', 'N', α, G, one(T), Φ)
        println("Block: $i, indices: $(indices[i])")
    end
    #Symmetric(Φ, :U) # Dette er et view så usikker.... Kan nok ikke subsettes! jo
    LinearAlgebra.copytri!(Φ, 'U')
end

# Load data
println("In julia")
gfile = ARGS[1]
sz = ARGS[2]
println("Input file: $(gfile)")
genedat_moba = SnpData(gfile)
mobadat = DataFrame(CSV.File(gfile * "_" * sz * ".moba", missingstring = "NA"))

# Step 1
# Filter plink data according to moba data
# --------------------------------------
# This is slow, but does not require much memory (~30 min first version genotype data)
# If we dont use des it will make .filtered fille in original storage location (genedat.src)
println("Filter genotype data according to moba data: ")
@time genedat_moba = SnpArrays.filter(genedat_moba, des = gfile * ".filtered", f_person = x -> x[:iid] in mobadat[!, :iid])
# Dest file might not work on cluster?
# Step 2
# Rearrange filtered plink data according to moba data
# --------------------------------------
println("Arrange genotype data according to moba data: ")
ord = something.(indexin(mobadat[!,:iid], genedat_moba.person_info[!,:iid])) # find index of filtered snp in moba
# This allocates
# 6 min all data
# I think you need to copt the file and open with r+
#cp(gfile * ".filtered.bim", gfile * ".filtered.reordered.bim", force = true)
#cp(gfile * ".filtered.bed", gfile * ".filtered.reordered.bed", force = true)
#cp(gfile * ".filtered.fam", gfile * ".filtered.reordered.fam", force = true)
#genedat_moba = SnpData(gfile * ".filtered.reordered", "r+")
@time SnpArrays.reorder!(genedat_moba, ord) # reorder snps
#@time SnpArrays.reorder!(a, ord) # reorder snps
genedat_moba.person_info
mobadat

# Step 3
# Compute grm on the filtered and arranged plink genotype data
# --------------------------------------
println("Compute grm: ")
# prøv dette
GC.gc()
@time A = 2 * grm_blocks(genedat_moba.snparray, 400)
write_grm(A, "results/grm$(size(A, 2))")