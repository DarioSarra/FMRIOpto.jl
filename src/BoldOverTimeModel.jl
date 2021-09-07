include("prepareDB.jl")
##
fdata = filter(r ->
    r.Area != "PIR" &&
    r.Area != "OT" &&
    r.Area != "ENT" &&
    r.Area != "EP" &&
    r.Area != "ENTm" &&
    r.Area != "ENTl" &&
    # r.Mouse != "FM11" &&
    # r.Mouse != "FM14" &&
    r.Valid_frs,
    # r.Area != "Ent" &&
    fulldata)
    # lvA = union(fdata.Area)
    # fdata.Area = categorical(fdata.Area; levels = lvA)
##
bold_df, plt, plt_df, to_keep, m1, mp, to_keep2, stim_df, m2, m3, p_simple, pcomplex, res =
    analysis_bold(fdata,:BOLD_85;
    method = :None, av_pres = true,
    peak_to = :all, nmax = 6, window = 5:8,
    nanmath = true, gamma = false,
    fconc = false, lim = (-0.4,1))
show(res)

#=
bold type
fconc = true use 3 highest concentration

method
    :Peak = from max onwards until peak_to & using the average response if av_pres = true
    :Max = which takes x nmax values from the average response per mouse and concentration
    :Window = average per animal of values in a window of frames (column Time)
=#
##
testdf = filter(r -> r.Area in to_keep2 && r.Area in to_keep && r.Genotype == "sert", bold_df)
# testdf = filter(r -> r.Area != "pPIR", bold_df)
# is Sertonin having an additive effect on bold?
f0 = @formula(Bold ~ 1 + Time + Concentration*Area + (1+Time|Mouse) + (1+Time|Area))
m0 = fit(MixedModel, f0,testdf)
f1 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim + (1+Time|Mouse) + (1+Time|Area))
m1 = fit(MixedModel, f1,testdf)
t1 = MixedModels.likelihoodratiotest(m0,m1)

# is the effect of Serotonin additive effect dependent on the area? No
f2 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim + Stim & Area + (1+Time|Mouse) + (1+Time|Area))
m2 = fit(MixedModel, f2,testdf)
t2 = MixedModels.likelihoodratiotest(m1,m2)

# is Serotonin having a multiplicative effect? YES
f3 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim + Stim & Concentration + (1+Time|Mouse) + (1+Time|Area))
m3 = fit(MixedModel, f3,testdf)
t3 = MixedModels.likelihoodratiotest(m1,m3)

# is the effect of Serotonin divisive effect dependent on the area? Yes
f4 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim + Stim & Concentration + Stim & Concentration & Area + (1+Time|Mouse) + (1+Time|Area))
m4 = fit(MixedModel, f4,testdf)
t4 = MixedModels.likelihoodratiotest(m3,m4)

# is a model wiht all interactions better?
f5 = @formula(Bold ~ 1 + Time + Concentration*Area*Stim + (1+Time|Mouse) + (1+Time|Area))
m5 = fit(MixedModel, f5,testdf)
t5 = MixedModels.likelihoodratiotest(m4,m5)

# is the divisive effect solely depending on area?
f6 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim + Stim & Concentration & Area + (1+Time|Mouse) + (1+Time|Area))
m6 = fit(MixedModel, f6,testdf)
t6 = MixedModels.likelihoodratiotest(m6,m4)

##
testdf = filter(r -> r.Genotype == "sert", bold_df)
# testdf = filter(r -> r.Area != "pPIR", bold_df)
# is Sertonin having an additive effect on bold?
f0 = @formula(Bold ~ 1 + Time + Concentration*Area + (1+Time|Mouse) + (1+Time|Area))
m0 = fit(MixedModel, f0,testdf)
f1 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim + (1+Time|Mouse) + (1+Time|Area))
m1 = fit(MixedModel, f1,testdf)
t1 = MixedModels.likelihoodratiotest(m0,m1)

# is Serotonin having a multiplicative effect? YES
f2 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim + Stim & Concentration + (1+Time|Mouse) + (1+Time|Area))
m2 = fit(MixedModel, f2,testdf)
t2 = MixedModels.likelihoodratiotest(m1,m2)

# is the effect of Serotonin additive effect dependent on the area? No
f3 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim  + Stim & Concentration + Stim & Area + (1+Time|Mouse) + (1+Time|Area))
m3 = fit(MixedModel, f3,testdf)
t3 = MixedModels.likelihoodratiotest(m2,m3)

# is the effect of Serotonin divisive effect dependent on the area? Yes
f4 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim + Stim & Concentration + Stim & Concentration & Area + (1+Time|Mouse) + (1+Time|Area))
m4 = fit(MixedModel, f4,testdf)
t4 = MixedModels.likelihoodratiotest(m2,m4)

# is a model wiht all interactions better?
f5 = @formula(Bold ~ 1 + Time + Concentration*Area*Stim + (1+Time|Mouse) + (1+Time|Area))
m5 = fit(MixedModel, f5,testdf)
t5 = MixedModels.likelihoodratiotest(m3,m5)
t5 = MixedModels.likelihoodratiotest(m4,m5)
##
testdf = filter(r -> r.Genotype == "sert", bold_df)
f5 = @formula(Bold ~ 1 + Time + Concentration*Area*Stim + (1+Time|Mouse) + (1+Time|Area))
m5 = fit(MixedModel, f5,testdf)
##
OFCdf = filter(r-> r.Genotype == "sert" && r.Area == "OFC", bold_df)
OBdf = filter(r-> r.Genotype == "sert" && r.Area == "OB", bold_df)
farea = @formula(Bold ~ 1 + Time + Concentration*Stim + (1+Time|Mouse))
mOFC = fit(MixedModel,farea,OFCdf)
mOB = fit(MixedModel,farea,OBdf)
