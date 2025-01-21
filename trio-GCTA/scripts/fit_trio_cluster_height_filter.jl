#author: Espen Eilertsen, Laura Hegemann 
#filter GRM on relatedness and fit full, no covariance, direct only, and null models for height phenotype


using VCModels
using SnpArrays, DataFrames, CSV, LinearAlgebra
using StatsBase, StatsModels
using JLD, HDF5

cd("/ess/p471/data/durable/people/Laura_H/Cluster_backup/trio_gcta_qc/")
#cd("Z:/projects/trio_gcta_qc")

function vec2ugrm(v::AbstractVector{T}) where {T <: AbstractFloat}
    k = length(v) # unique elements / triangle + diagonal
    p = Int((sqrt(8k + 1) - 1) / 2) # rows/cols in matrix
    G = zeros(T, p, p)
    l = 1
    @inbounds for i in 1:p # cols
        @inbounds for j in 1:i # rows
            G[j, i] = v[l]
            l += 1
        end
    end
    G
end

function read_grm(filename::AbstractString) 
    nme = "grm_u"
    g = h5open(filename * ".h5", "r") do file
        read(file, nme)
    end
    # Symmetric(vec2ugrm(g), :U)
    LinearAlgebra.copytri!(vec2ugrm(g), 'U')
    #vec2ugrm(g)
end

function findGRMs(A, c_inds_use, m_inds_use, f_inds_use) 
    # Grms
    A_m = A[m_inds_use, m_inds_use]
    A_f = A[f_inds_use, f_inds_use]
    A_c = A[c_inds_use, c_inds_use]
    D_fm = A[f_inds_use, m_inds_use] + A[m_inds_use, f_inds_use]
    D_cm = A[c_inds_use, m_inds_use] + A[m_inds_use, c_inds_use]
    D_cf = A[c_inds_use, f_inds_use] + A[f_inds_use, c_inds_use]
    R = Diagonal(ones(size(A_c, 1)))
    [Symmetric(A_m), Symmetric(A_f), Symmetric(A_c), Symmetric(D_fm), Symmetric(D_cm), Symmetric(D_cf), R]
end

function sel_trio(A, m_inds, f_inds, c_inds)
     # Check for relatives
    mf_inds = vcat(m_inds, f_inds)
    A_mf = A[mf_inds, mf_inds]
     # A_mf = view(A, mf_inds, mf_inds) # This view may be wrong and remove parent-offspring!
    keep_A_mf = kinship_pruning(A_mf, method=:bottom_up, cutoff=0.10)
    m_inds_keep = keep_A_mf[m_inds]
    f_inds_keep = keep_A_mf[f_inds]
    mf_inds_keep = findall((m_inds_keep + f_inds_keep) .== 2)
     
    m_inds_use = m_inds[mf_inds_keep]
    f_inds_use = f_inds[mf_inds_keep]
    c_inds_use = c_inds[mf_inds_keep]
    (m_inds_use, f_inds_use, c_inds_use)
end



# Load data
sz = "107088"
gfile = "MoBaPsychGen_v1_filter_info_0.95" * "_" * sz
mobadat = DataFrame(CSV.File(gfile * ".moba", missingstring="NA"))
pcdat = DataFrame(CSV.File("/ess/p471/cluster/data/genetic_data/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt", missingstring="NA"))
A = read_grm("grm107088")

mobadat.order = 1:size(mobadat, 1)
mobadat = leftjoin(mobadat, pcdat, on = [:iid => :IID])
sort!(mobadat, :order)
# Subset grm blocks, after possibly filtering high relatedness coefficients
# --------------------------------------
# This is fast
A = Symmetric(A, :U)

m_inds = 1:35696
f_inds = 35697:71392
c_inds = 71393:107088

# Check pos
ii = 14000
trio_inds_i = vcat(m_inds[ii], f_inds[ii], c_inds[ii])
println(A[trio_inds_i, trio_inds_i])

m_inds_sel, f_inds_sel, c_inds_sel = sel_trio(A, m_inds, f_inds, c_inds)
trio_inds_sel_i = vcat(m_inds_sel[ii], f_inds_sel[ii], c_inds_sel[ii])
println(A[trio_inds_sel_i, trio_inds_sel_i])

# Step 5
# Fit models
# --------------------------------------
mobadat[!, :sex] .= string.(mobadat[!, :sex])
mobadat[!, :Age_height] .= string.(mobadat[!, :Age_height])
mobadat[!, :Age_exam] .= string.(mobadat[!, :Age_exam])
mobadat[!, :genotyping_batch_num] .= string.(mobadat[!, :genotyping_batch_num])


vari = "height"

pos_use = findall(!ismissing, mobadat[c_inds_sel, vari]) 
println(length(pos_use))
mobadat[c_inds_sel[pos_use],vari]

m_inds_sel_nomiss = m_inds_sel[pos_use]
f_inds_sel_nomiss = f_inds_sel[pos_use]
c_inds_sel_nomiss = c_inds_sel[pos_use]

m_pc = Matrix(pcdat[m_inds_sel_nomiss, Symbol.(:PC, 1:20)])
f_pc = Matrix(pcdat[f_inds_sel_nomiss, Symbol.(:PC, 1:20)])
mf_pc = m_pc + f_pc
    #y = convert.(Float64, mobadat[c_inds_sel_nomiss, vari])
mf = ModelFrame(term(Symbol(vari)) ~ ConstantTerm(1) + term(:genotyping_batch_num) + term(:sex) + term(:Age_height), mobadat[c_inds_sel_nomiss, :])
X = ModelMatrix(mf).m
X = hcat(X, mf_pc)
y = response(mf)
y = float(y)

rs = findGRMs(A, c_inds_sel_nomiss, m_inds_sel_nomiss, f_inds_sel_nomiss)


# Full model
# -------------------------------------------
    function VCModels.transform!(δ::Vector, θ::Vector)
        δ[1] = θ[1]^2    
        δ[2] = θ[2]^2 + θ[4]^2    
        δ[3] = θ[3]^2 + θ[5]^2 + θ[6]^2    
        δ[4] = θ[2] * θ[1]
        δ[5] = θ[3] * θ[1]
        δ[6] = θ[3] * θ[2] + θ[5] * θ[4]
        δ[7] = θ[7]
        println("transformation performed")
    end
    dat = VCData(y, X, rs)
    lbs = [0.0, -Inf, -Inf, 0.0, -Inf, 0.0, 0.0] # not sure what this does - does it need to be changed if focal individual changes?
    vd = var(y)
    
    ini = sqrt.(vd .* [0.1, 0.0, 0.0, 0.1, 0.0, 0.1, 0.7]) #are these starting params? changed since child is focal ind. 
    @time mod = VCModel(dat, ini, lbs, false)
    mod.opt.ftol_abs = 0
    mod.opt.ftol_rel = 0 
    @time VCModels.fit!(mod)
    hessian!(mod.opt.H, mod)
    JLD.save("results/" * vari * "_full_filter_update_$(size(mobadat, 1)).jld", "model", mod)

# No covariances
# -------------------------------------------
    function VCModels.transform!(δ::Vector, θ::Vector)
        δ[1] = θ[1]^2    
        δ[2] = θ[2]^2    
        δ[3] = θ[3]^2  
        δ[4] = θ[4]
        println("transformation performed")
    end
    dat = VCData(y, X, [rs[1], rs[2], rs[3], rs[7]])
    lbs = [-Inf, -Inf, 0.0, 0.0] #wasn't sure about this #could be 0
    vd = var(y)
    ini = sqrt.(vd .* [0.0, 0.0, 0.3, 0.7])
    @time mod = VCModel(dat, ini, lbs, false)
    @time VCModels.fit!(mod)
    hessian!(mod.opt.H, mod)
    JLD.save("results/" * vari * "_nocov_filter_update_$(size(mobadat, 1)).jld", "model", mod)

# direct effects only
function VCModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]
    println("transformation performed")
end
dat = VCData(y, X, [rs[3], rs[7]])
lbs = [0.0, 0.0]
ini = sqrt.(vd .* [0.10, 0.9^2])

@time mod = VCModel(dat, ini, lbs, false)
@time fit!(mod)
hessian!(mod.opt.H, mod)
JLD.save("results/" * vari * "_direct_filter_update_$(size(mobadat, 1)).jld", "model", mod)


#Null model 

function VCModels.transform!(δ::Vector, θ::Vector)
δ[1] = θ[1]
end
dat = VCData(y, X, [rs[7]])
lbs = [0.0]
ini = sqrt.(vd .* [1.0])

@time mod = VCModel(dat, ini, lbs, false)
@time fit!(mod)
hessian!(mod.opt.H, mod)
JLD.save("results/" * vari * "_null_filter_update_$(size(mobadat, 1)).jld", "model", mod)
