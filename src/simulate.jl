"""
    function simu(; nid = 100, nlc = 1000, nql = 100)
---
## Simulate:
- 1000 loci
- 100  ID
- 100  QTL
- Breeding values
- VR-I G
"""
function simu()
    gtxt = joinpath(ddat, "gt.txt")
    nid, nlc, nql = 100, 1000, 100
    Z   = convert.(Float64, reshape(read(gtxt), nlc, :)) .- 48
    qlc = sort(rand(1:nlc, nql))       # sample QTL
    eql = randn(nql) ./ sqrt(nql)      # QTL effects
    ID  = String[]                     # random ID names
    for _ in 1:nid
        push!(ID, randstring(6))
    end
    qtl = Z[qlc, :]
    bv = qtl'eql                # true breeding values
    twop = mean(Z, dims=2)      # for G matrix
    Z  .-= twop
    s2pq = (1 .- .5twop)'twop
    r2pq = 1 / s2pq[1]          # mul faster than div
    G    = Z'Z .* r2pq + 0.0001I
    
    ID, bv, G
end    

"""
    function homo_R(bv; h2 = 0.8)
---
Simulate phenotypes with homogeneous residuals.
"""
function homo_R(bv; h2 = 0.8)
    nid = length(bv)
    vg  = var(bv)
    ve  = vg/h2*(1-h2)
    se  = sqrt(ve)
    pp  = bv + randn(nid) .* se
end

"""
    function hete_R(bv; h2 = 0.8, mo = 20)
---
Suppose ID may have 1:mo records
"""
function hete_R(bv; h2 = 0.8, mo = 20)
    nid = length(bv)
    nob = rand(1:mo, nid)
    vg  = var(bv)
    ve  = vg/h2*(1-h2)
    se  = sqrt(ve)
    wt  = 1 ./ sqrt.(nob)
    pp  = bv + randn(nid) .* wt
    nob, pp
end
