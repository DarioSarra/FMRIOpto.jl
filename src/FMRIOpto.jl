module FMRIOpto

    using Reexport
    @reexport using DataFrames, CSV, CategoricalArrays, BrowseTables
    @reexport using Statistics, StatsBase, Plots, StatsPlots, HypothesisTests
    @reexport using MixedModels

include("functions.jl")

export nanmean, nansem, safemean, safesem
export plot_bold, defineBold

end # module
