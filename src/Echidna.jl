module Echidna

using Serialization, Statistics, Random, LinearAlgebra

dpkg = joinpath(splitpath(pathof(Echidna))[1:end-2]...)
ddat = joinpath(dpkg, "dat")
dtmp = joinpath(dpkg, "tmp")
isdir(dtmp) || mkdir(dtmp)

include("workflow.jl")
include("simulate.jl")
include("compare.jl")

export workflow

end # module
