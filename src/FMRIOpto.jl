module FMRIOpto

    using Reexport
    @reexport using DataFrames, CSV, CategoricalArrays, BrowseTables
    @reexport using Statistics, StatsBase, Plots, StatsPlots, HypothesisTests
    @reexport using MixedModels

include("functions.jl")

export analyze_bold, max_bold, average_boldtime

end # module
