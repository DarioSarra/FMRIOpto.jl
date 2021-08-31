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
