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
sort(fulldata.Area)
lvA = union(fulldata.Stim)
# fulldata.Area = categorical(fulldata.Area; levels = lvA)
# lvC = union(fulldata.Concentration)
# fulldata.Concentration = categorical(fulldata.Concentration; levels = lvC)
########
des = describe(fulldata)
open_html_table(des)
des2 = combine(groupby(fulldata, :Mouse),
    [:Genotype, :Area, :Concentration, :Stim, :Repetition_nr, :Trial_nr, :Block] .=> (c-> [union(c)]), nrow)
open_html_table(des2)
## average response per concentraition and stim over all presentation
fdata = filter(r ->
    r.Area != "PIR" &&
    # r.Area != "Ent" &&
    r.Valid_frs,
    # r.Mouse != "FM14",
    fulldata)
gd0 = groupby(fdata,[:Genotype, :Mouse, :Stim,:Area, :Concentration, :Time])
df0 = combine(gd0, [:BOLD_90, :BOLD_85,:BOLD_0] .=> mean .=> [:BOLD_90, :BOLD_85,:BOLD_0])
##
fdata2 = filter(r -> r.Genotype == "sert" &&
    r.Concentration != 0.0001,
    # r.Mouse != "FM14",
    df0)

gd1 = groupby(fdata2, [:Mouse, :Area, :Stim, :Concentration])



df1 = combine(gd1, [:BOLD_90, :BOLD_85,:BOLD_0] .=> max_bold .=> [:BOLD90_max6, :BOLD85_max6,:BOLD0_max6],
    [:BOLD_90, :Time] => ((b,t) -> average_boldtime(b,t)) => :BOLD90_4to10,
    [:BOLD_85, :Time] => ((b,t) -> average_boldtime(b,t)) => :BOLD85_4to10,
    [:BOLD_0, :Time] => ((b,t) -> average_boldtime(b,t)) => :BOLD0_4to10,
    )
## Areas that respond to concentration
df2 = filter(r -> !r.Stim, df1)
vaform2c = @formula(BOLD90_4to10 ~ 1 + Concentration + (1|Mouse))
mm2c = fit(MixedModel, vaform2c,df2)
vaform2 = @formula(BOLD90_4to10 ~ 1 + Concentration*Area + (1|Mouse))
mm2 = fit(MixedModel, vaform2,df2)
MixedModels.likelihoodratiotest(mm2c,mm2)

res = DataFrame(Area = coefnames(mm2), P = mm2.pvalues, Z = mm2.beta ./ mm2.stderror)
res2 = filter!(r -> occursin("Concentration & Area: ", r.Area), res)
res2.Area = replace.(res2.Area, "Concentration & Area: " => "")
sort!(res2,:P)
res3 = filter(r -> r.P <= 0.049, res2)
##
df3 = filter(r -> r.Area in res3.Area, df1)
vaform3c = @formula(BOLD90_4to10 ~ 1 + Concentration*Stim + Area*Stim + (1|Mouse))
mm3c = fit(MixedModel, vaform3c,df3)
vaform3 = @formula(BOLD90_4to10 ~ 1 + Concentration*Area*Stim + (1|Mouse))
mm3 = fit(MixedModel, vaform3,df3)
MixedModels.likelihoodratiotest(mm3c,mm3)


# show(MIME("text/latex"), mm1)
##
Area = "OB"
df4 = filter(r -> r.Area == Area, df1)
vaform4c = @formula(BOLD90_4to10 ~ 1 + Concentration + Stim + (1|Mouse))
mm4c = fit(MixedModel, vaform4c,df4)
vaform4 = @formula(BOLD90_4to10 ~ 1 + Concentration*Stim + (1|Mouse))
mm4 = fit(MixedModel, vaform4,df4)
MixedModels.likelihoodratiotest(mm4c,mm4)
##
Area = "OFC"
df5 = filter(r -> r.Area == Area, df1)
vaform5c = @formula(BOLD85_max6 ~ 1 + Concentration + Stim + (1|Mouse))
mm5c = fit(MixedModel, vaform5c,df5)
vaform5 = @formula(BOLD85_max6 ~ 1 + Concentration*Stim + (1|Mouse))
mm5 = fit(MixedModel, vaform5,df5)
MixedModels.likelihoodratiotest(mm5c,mm5)
##
fdata = filter(r ->
    r.Area != "PIR" &&
    r.Area != "OT" &&
    # r.Area != "Ent" &&
    r.Valid_frs,
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
##
Area = "OB"
df4 = filter(r -> r.Area == Area, df1)
form_4c = @formula(BOLD_90 ~ 1 + Time + Concentration + Stim + (1+Time|Mouse))
model_4c = fit(MixedModel, form_4c,df4)
form_4 = @formula(BOLD_90 ~ 1 + Time + Concentration*Stim + (1+Time|Mouse))
model_4 = fit(MixedModel, form_4,df4)
MixedModels.likelihoodratiotest(model_4c,model_4)
##
Area = "OFC"
df5 = filter(r -> r.Area == Area, df1)
form_5c = @formula(BOLD_90 ~ 1 + Time + Concentration + Stim + (1+Time|Mouse))
model_5c = fit(MixedModel, form_5c,df5)
form_5 = @formula(BOLD_90 ~ 1 + Time + Concentration*Stim + (1+Time|Mouse))
model_5 = fit(MixedModel, form_5,df5)
MixedModels.likelihoodratiotest(model_5c,model_5)
##
Area = "aPIR"
df5 = filter(r -> r.Area == Area, df1)
form_5c = @formula(BOLD_90 ~ 1 + Time + Concentration + Stim + (1+Time|Mouse))
model_5c = fit(MixedModel, form_5c,df5)
form_5 = @formula(BOLD_90 ~ 1 + Time + Concentration*Stim + (1+Time|Mouse))
model_5 = fit(MixedModel, form_5,df5)
MixedModels.likelihoodratiotest(model_5c,model_5)

model_5


##
vaform2 = @formula(BOLD_90 ~ 1 + Concentration*Area + (1|Mouse));
df2 = filter(r -> !r.Stim, df1)
mm2 = fit(MixedModel, vaform2,df2)

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
