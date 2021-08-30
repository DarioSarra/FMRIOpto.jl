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
    :BOLD_0 => Float64
    )
using FMRIOpto
fulldata =  CSV.read(joinpath(files_loc,"odour_dose_data.csv"), DataFrame;
    types = columns_types,
    missingstring = "NaN",
    limit = 14080)
lvA = union(fulldata.Area)
fulldata.Area = categorical(fulldata.Area; levels = lvA)
# lvC = union(fulldata.Concentration)
# fulldata.Concentration = categorical(fulldata.Concentration; levels = lvC)
########
des = describe(fulldata)
open_html_table(des)
des2 = combine(groupby(fulldata, :Mouse),
    [:Genotype, :Area, :Concentration, :Stim, :Repetition_nr, :Trial_nr, :Block] .=> (c-> [union(c)]), nrow)
open_html_table(des2)
##
vaform = @formula(BOLD_90 ~ 1 + Concentration*Area*Stim + (1|Mouse+Area));
fdata = filter(r -> r.Genotype == "sert" &&
    r.Area != "PIR" &&
    r.Concentration != 0.0001 &&
    r.Mouse != "FM14",
    fulldata)
mm1 = fit(MixedModel, vaform,fdata)
show(MIME("text/latex"), mm1)
open_html_table(filter(r -> r.Mouse == "FM12", fulldata))
##
vaform2 = @formula(BOLD_90 ~ 1 + Concentration*Stim + (1|Mouse));
cnames = ["z(Intercept)",
    "zConcentration",
    "zStim",
    "zConcentration & Stim",
    "p(Intercept)",
    "pConcentration",
    "pStim",
    "pConcentration & Stim"]
allmodels = combine(groupby(fdata, :Area)) do dd
    try
        mm = fit(MixedModel, vaform2, dd)
        DataFrame(
            vcat("z".*coefnames(mm),"p".*coefnames(mm)) .=>
            vcat(mm.beta ./ mm.stderror,
            mm.pvalues)
            )
    catch
        # (mm = missing,)
        DataFrame(
            cnames .=>
            missing
            )
    end
end
open_html_table(allmodels)
##
coef(allmodels.mm[6])
fixef(allmodels.mm[6])
propertynames(allmodels.mm[6])

coefnames(allmodels.mm[6])
allmodels.mm[6].beta
allmodels.mm[6].stderror
allmodels.mm[6].beta ./ allmodels.mm[6].stderror
allmodels.mm[6].pvalues
DataFrame(
    vcat("z".*coefnames(allmodels.mm[6]),"p".*coefnames(allmodels.mm[6])) .=>
    vcat(allmodels.mm[6].beta ./ allmodels.mm[6].stderror,
    allmodels.mm[6].pvalues)
    )
DataFrame(
    vcat("z".*coefnames(allmodels.mm[6]),"p".*coefnames(allmodels.mm[6])) .=>
    missing
    )

vcat("z".*coefnames(allmodels.mm[6]),"p".*coefnames(allmodels.mm[6]))
