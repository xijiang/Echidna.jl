"""
    function comp_homo()
---
Compare EBV with homogeneous residuals of MME and Echidna results.
I use the last 10 ID as validators.
"""
function comp_homo()
    ID, bv, G = simu()
    pp = homo_R(bv)
    ebv = 1
end
