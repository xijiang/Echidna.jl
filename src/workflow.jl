function workflow(; nv = 10)
    ########################################
    println("Simulation")
    println("  - ID, breeding values, and G-matrix")
    ID, bv, G = simu()
    giv = inv(G)
    println("  - homogeneous R")
    ya = homo_R(bv)
    println("  - heterogeneous R")
    nob, yb = hete_R(bv)

    ########################################
    println("\nPrediction with homo data")
    # model: y ~ mu + ebv + e
    # using the last 10 ID for validation
    nid, nv = length(ID), 10
    nt = nid - nv
    println("  - Prepare matrices")
    X = ones(nt)
    Z = [I zeros(nt, nv)]
    println("  - for the homogeneous R data")
    λ = 0.25                    # (1-h2)/h2
    lhs = [ nid X'Z
            Z'X (Z'Z + giv .* λ) ]
    y = ya[1:nt]
    rhs = [ X'y
            Z'y ]
    ebv = lhs \ rhs
    cor_1 = cor(bv[nt+1:end], ebv[nt+2:end])

    println("\nUsing Echidna")
    
end
