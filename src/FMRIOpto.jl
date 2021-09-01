module FMRIOpto

    using Reexport
    @reexport using DataFrames, CSV, CategoricalArrays, BrowseTables
    @reexport using Statistics, StatsBase, Plots, StatsPlots, HypothesisTests
    @reexport using MixedModels

include("utilities.jl")
include("functions.jl")

export nanmean, nansem, safemean, safesem, iscolumn, model_res
export average_presentation, define_bold, plot_bold, odour_response, stim_interaction

end # module
