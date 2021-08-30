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
##
fdata = filter(r ->
    r.Area != "PIR" &&
    r.Area != "OT" &&
    r.Valid_frs,
    # r.Area != "Ent" &&
    # r.Mouse != "FM14",
    fulldata)
    # lvA = union(fdata.Area)
    # fdata.Area = categorical(fdata.Area; levels = lvA)
gd0 = groupby(fdata,[:Genotype, :Mouse, :Stim, :Area, :Concentration])
df0 = combine(gd0) do dd
    any(ismissing,dd.BOLD_0) && println("got NaN ", dd.Mouse[1], dd.Area[1], dd.Stim[1], dd.Concentration[1])
    if all(ismissing.(dd.BOLD_0))
        return filter(r -> r.Area == "NoData",dd)
    else
        i = findmax(skipmissing(dd.BOLD_0))[2]
        return dd[i:end,:]
    end
end
##
f0 = @formula(Bold ~ 1 + Time + Concentration + Area + (1+Time|Mouse+Area))
f1 = @formula(Bold ~ 1 + Time + Concentration + Area + Stim + Concentration & Area + Concentration & Area & Stim + (1+Time|Mouse))
f2 = @formula(Bold ~ 1 + Time + Concentration*Area*Stim + (1+Time|Mouse+Area))
df1 = filter(r -> !(r.Area in ["pPIR", "OT", "NDB", "LHA", "COA"]) &&
    r.Concentration != 0.0001,
    # r.Mouse != "FM14",
    df0)
df1.Bold = df1.BOLD_90
m0 = fit(MixedModel, f0,df1)
m1 = fit(MixedModel, f1,df1)
m2 = fit(MixedModel, f2,df1)
##
df1 = filter(r -> r.Genotype == "sert" &&
    r.Concentration != 0.0001,
    # r.Mouse != "FM14",
    df0)
form_1c = @formula(BOLD_0 ~ 1 + Time + Concentration*Area + Stim*Area + (1+Time|Mouse+Area))
model_1c = fit(MixedModel, form_1c,df1)
form_1 = @formula(BOLD_0 ~ 1 + Time + Concentration*Area*Stim + (1+Time|Mouse+Area))
model_1 = fit(MixedModel, form_1,df1)
MixedModels.likelihoodratiotest(model_1c,model_1)
## Areas that respond to concentration
df2 = filter(r -> !r.Stim, df0)
form_2c = @formula(BOLD_90 ~ 1 + Time + Concentration + (1+Time|Mouse))
model_2c = fit(MixedModel, form_2c,df2)
form_2 = @formula(BOLD_90 ~ 1 + Time + Concentration*Area + (1+Time|Mouse))
model_2 = fit(MixedModel, form_2,df2)
MixedModels.likelihoodratiotest(model_2c,model_2)
res = DataFrame(Area = coefnames(model_2), P = model_2.pvalues, Z = model_2.beta ./ model_2.stderror)
res2 = filter!(r -> occursin("Concentration & Area: ", r.Area), res)
res2.Area = replace.(res2.Area, "Concentration & Area: " => "")
sort!(res2,:P)
res3 = filter(r -> r.P <= 0.049, res2)
##
df3 = filter(r -> r.Area in res3.Area, df1)
union(df3.Area)
df3.Area = categorical(df3.Area; levels = union(df3.Area))

form_3c = @formula(BOLD_90 ~ 1 + Time + Concentration*Stim + Area*Stim + (1+Time|Mouse))
model_3c = fit(MixedModel, form_3c,df3)

form_3 = @formula(BOLD_90 ~ 1 + Time + Concentration*Area*Stim + (1+Time|Mouse))
model_3 = fit(MixedModel, form_3,df3)

form_3b = @formula(BOLD_90 ~ 1 + Time + Concentration + Area + Stim + Concentration & Area + Concentration & Area & Stim + (1+Time|Mouse))
model_3b = fit(MixedModel, form_3b,df3)
MixedModels.likelihoodratiotest(model_3b,model_3)
