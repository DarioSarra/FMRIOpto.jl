include("prepareDB.jl")
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
##
df0 = defineBold(fdata, :BOLD_0; method = :Window, pick_to = :all, nmax = 4, window = 4:10, nanmath = true)

## plot areas
plt, meandf = plot_bold(df0, :Bold; nanmath = true)
plot(plt..., size = (900*1.5,1200*1.5))#layout = (6,3)
## Areas that respond to concentration
df1 = filter(r -> !r.Stim &&
    !(r.Area in ["pPIR", "COA","LHA"]), df0)
f0 = @formula(Bold ~ 1 + Time + Concentration + (1+Time|Mouse))
f1 = @formula(Bold ~ 1 + Time + Concentration*Area + (1+Time|Mouse))
df1.Bold = df1.BOLD_90
m0 = fit(MixedModel, f0, df1)
m1 = fit(MixedModel, f1, df1)
MixedModels.likelihoodratiotest(m0,m1)
res = DataFrame(Area = coefnames(m1), P = m1.pvalues, Z = m1.beta ./ m1.stderror)
res2 = filter!(r -> occursin("Concentration & Area: ", r.Area), res)
res2.Area = replace.(res2.Area, "Concentration & Area: " => "")
sort!(res2,:P)
res3 = filter(r -> r.P <= 0.049, res2)
##
f3 = @formula(Bold ~ 1 + Time + Concentration*Area + (1+Time|Mouse))
f4 = @formula(Bold ~ 1 + Time + Concentration + Area + Stim + Concentration & Area + Concentration & Area & Stim + (1+Time|Mouse))
f5 = @formula(Bold ~ 1 + Time + Concentration*Area*Stim + (1+Time|Mouse))
df2 = filter(r -> r.Area in res3.Area &&
    r.Concentration != 0.0001,
    # r.Mouse != "FM14",
    df0)
m3 = fit(MixedModel, f3,df2)
m4 = fit(MixedModel, f4,df2)
m5 = fit(MixedModel, f5,df2)
MixedModels.likelihoodratiotest(m3,m4)
MixedModels.likelihoodratiotest(m4,m5)
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
