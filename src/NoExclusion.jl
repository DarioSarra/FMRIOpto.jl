const mac_gdrive = "/Volumes/GoogleDrive/My Drive"
const linux_gdrive = "/home/beatriz/mainen.flipping.5ht@gmail.com"
const files_dir = "Flipping/Datasets/FMRIOpto"
const figs_dir = "Flipping/Datasets/FMRIOpto/figures/"
## load
if ispath(linux_gdrive)
    ongoing_dir = linux_gdrive
else
    ongoing_dir = mac_gdrive
end

files_loc = joinpath(ongoing_dir,files_dir)
figs_loc = joinpath(ongoing_dir,figs_dir)

const columns_types = Dict(
    :Mouse => String,
    :Genotype => String,
    :Area => String,
    :Concentration => Float64,
    :Stim => Bool,
    :Repetition_nr => Int,
    :Trial_nr => Int,
    :Block => Int,
    :BOLD_90 => Float64,
    :BOLD_85 => Float64,
    :BOLD_0 => Float64,
    :Frame => Int,
    :Valid_frs => Bool,
    :Outlier_count_per_rep => Int
    )
using FMRIOpto
fulldata =  CSV.read(joinpath(files_loc,"odour_dose_data.csv"), DataFrame;
    types = columns_types,
    missingstring = "NaN")
mapcols!(Array, fulldata)
fulldata.Time = fulldata.Frame .+ 3
## average response per concentraition and stim over all presentation
gd0 = groupby(fulldata,[:Genotype, :Mouse, :Stim,:Area, :Concentration, :Time])
df0 = combine(gd0, [:BOLD_90, :BOLD_85,:BOLD_0] .=> mean .=> [:BOLD_90, :BOLD_85,:BOLD_0])

## define bold
fdata0 = filter(r -> r.Genotype == "sert" &&
    !(r.Area in ["pPIR", "OT", "NDB", "LHA", "COA"]) &&
    r.Concentration != 0.0001,
    # r.Mouse != "FM14",
    df0)
gd1 = groupby(fdata0, [:Mouse, :Area, :Stim, :Concentration])
df1 = combine(gd1, :BOLD_90 => mean => :BOLD_90)
plt, res, m, mOFC, mOB = analyze_bold(df1,:BOLD_90)
plt
res
m
mOFC
mOB
##
df2 = combine(groupby(fdata0, [:Genotype, :Mouse, :Stim,:Area, :Concentration]),
    [:BOLD_90, :BOLD_85,:BOLD_0] .=> max_bold .=> [:BOLD_90, :BOLD_85,:BOLD_0])
plt, res, m, mOFC, mOB, maPIR = analyze_bold(df2,:BOLD_90)
plt
res
m
mOFC
mOB
