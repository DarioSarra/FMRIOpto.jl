include("prepareDB.jl")
##
fdata = filter(r ->
    r.Area != "PIR" &&
    r.Area != "ENT" &&
    r.Area != "OT" &&
    r.Areac!=
    # r.Mouse != "FM14",
    r.Valid_frs,
    fulldata)
    # lvA = union(fdata.Area)
    # fdata.Area = categorical(fdata.Area; levels = lvA)
##
function prov_define_bold(df0, signal; av_pres = false, method = :Peak, peak_to = :all, nmax = 4, window = 4:10, nanmath = true)
    fmean, fsem = nanmath ? (safemean, safesem) : (mean, sem)
    if  method == :Max
        #average response over all presentation
        df1 = average_presentation(df0,signal; nanmath = nanmath)
        gd1 = groupby(df1,[:Genotype, :Mouse, :Stim, :Area, :Concentration])
        signal = :Bold
        #find the first n max values in the bold signal
        df2 = combine(gd1, signal => (s -> FMRIOpto.max_bold(s; max = nmax, nanmath = nanmath)) => :Bold)
    elseif method == :Window
        #average all bold values per mouse and concentration in a time window
        df1 = combine(groupby(df0, [:Time,:Concentration,:Area,:Stim,:Mouse,:Genotype]), signal => (s -> mean(skipmissing(s))) => :Bold)
        filter!(r -> r.Time in window, df1)
        gd1 = groupby(df1,[:Genotype, :Mouse, :Stim, :Area, :Concentration])
        df2 = combine(gd1, :Bold => (s -> mean(skipmissing(s))) => :Bold)
    elseif method == :None
        df2 = filter(r -> r.Time in window, df0)
        df2.Bold = df2[:,signal]
    else
        error("unrecognised bold processing method:/n options are :Peak, :Max, :Timewindow")
    end
    return df2
end

function prov_plot_bold(df0; nanmath = true, lim = (-0.4,0.7))
    fmean = nanmath ? safemean : mean
    fsem = nanmath ? safesem : sem
    #calculate mean repsonse over concentrations per area
    df = combine(groupby(df0,[:Area, :Stim, :Concentration, :Genotype]),  :Bold .=> [fmean, fsem] .=> [:mean, :sem])
    df1 = filter(r-> r.Genotype == "sert",df)
    # identify regions with missing concentration info
    res = combine(groupby(df1,[:Area]), :mean => (m -> .!(any(isnan.(m)))) => :to_keep)
    to_keep = filter(r -> r.to_keep, res)[:,:Area]
    # prepare a grid plot with the results for each regions
    df1.Color = [x ? :red : :black for x in df1.Stim]
    # posdict = Dict(v=>p for (v,p) in zip(union(df1.Concentration), collect(1:length(union(df1.Concentration)))))
    # df1.Pos = [get(posdict, x,0) for x in df1.Concentration]
    # sort!(df1,[:Stim,:Pos])
    sort!(df1,[:Stim,:Concentration])
    # fconc && filter!(r -> r.Concentration != 0.0001 ,df1)
    plt = []
        combine(groupby(df1, :Area)) do dd
            pp = @df dd scatter(:Concentration, :mean, yerror = :sem, group = :Stim, color = :Color,
                title = :Area[1], legend = false, ylims = lim)
            @df dd plot!(:Concentration, :mean, color = :Color, group = :Stim)
            Plots.abline!(0,0,color = :black, linestyle = :dash)
            push!(plt,pp)
        end
    return plt, df, to_keep
end
function prov_odour_response(df0; farea = :none)
    # use only control data to determine odour concentration dose dependency
    df1 = filter(r -> !r.Stim, df0)
    #remove lowest concentration
    # fconc && filter!(r -> r.Concentration != 0.0001, df1)
    # filters area with insufficient response, identified from plot_bold function
    farea != :none && filter!(r -> r.Area in farea, df1)
    # identify if the model requires to incorporate time
    if iscolumn(df0,:Time)
        f0 = @formula(Bold ~ 1 + Time + Concentration + (1+Time|Mouse))
        f1 = @formula(Bold ~ 1 + Time + Concentration*Area + (1+Time|Mouse))
    else
        f0 = @formula(Bold ~ 1 + Concentration + (1|Mouse))
        f1 = @formula(Bold ~ 1 + Concentration*Area + (1|Mouse))
    end
    # #replace NaN values used for plotting to missing for model
    # df1.Bold = [isnan(x) ? missing : x for x in df1.Bold]
    # model fit and test
    m0 = fit(MixedModel, f0, df1)
    m1 = fit(MixedModel, f1, df1)
    mp = MixedModels.likelihoodratiotest(m0,m1)
    # identifies area that significantly respond to odour concentration
    res = DataFrame(Area = coefnames(m1), P = m1.pvalues, Z = m1.beta ./ m1.stderror)
    res2 = filter!(r -> occursin("Concentration & Area: ", r.Area), res)
    res2.Area = replace.(res2.Area, "Concentration & Area: " => "")
    to_keep = filter(r -> r.P <= 0.049, res2)[:,:Area]

    return m1, mp, to_keep
end

function prov_stim_interaction(df0, farea)
    # filters area lacking concentration response identified by odour_response function
    df1 = filter(r -> r.Genotype == "sert" && r.Area in farea, df0)
    # fconc && filter!(r -> r.Concentration != 0.0001, df0)
    # identify if the model requires to incorporate time
    if iscolumn(df0,:Time)
        f1 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim + zerocorr(1+Time+Concentration|Mouse) + zerocorr(1+Time+Concentration|Area))
        f2 = @formula(Bold ~ 1 + Time + Concentration + Area + Stim + Concentration & Area + Concentration & Area & Stim + zerocorr(1+Time+Concentration|Mouse) + zerocorr(1+Time+Concentration|Area))
        f3 = @formula(Bold ~ 1 + Time + Concentration*Area*Stim + zerocorr(1+Time+Concentration|Mouse) + zerocorr(1+Time+Concentration|Area))
    else
        f1 = @formula(Bold ~ 1 + Concentration*Area + Stim + zerocorr(1+Concentration|Mouse) + zerocorr(1+Concentration|Area))
        f2 = @formula(Bold ~ 1 + Concentration + Area + Stim + Concentration & Area + Concentration & Area & Stim + zerocorr(1+Concentration|Mouse) + zerocorr(1+Concentration|Area))
        f3 = @formula(Bold ~ 1 + Concentration*Area*Stim + zerocorr(1+Concentration|Mouse) + zerocorr(1+Concentration|Area))
    end
    #replace NaN values used for plotting to missing for model
    # df1.Bold = [isnan(x) ? missing : x for x in df1.Bold]
    m1 = fit(MixedModel, f1,df1)
    m2 = fit(MixedModel, f2,df1)
    m3 = fit(MixedModel, f3,df1)
    p_simple = MixedModels.likelihoodratiotest(m1,m2)
    p_complex = MixedModels.likelihoodratiotest(m2,m3)
    res = model_res(m2)
    filter!(r -> occursin("Stim", r.Area), res)
    sort!(res,:P)
    return df1, m2, m3, p_complex, res
end

function prov_analysis_bold(df0, signal; fconc = true, av_pres = true, method = :Peak, peak_to = :all, nmax = 6, window = 4:12, logarithmic = false, nanmath = true, lim = (-0.4,0.7))
    #adjust bold signal
    bold_df = prov_define_bold(df0,signal; method = method, av_pres = av_pres, peak_to = peak_to, nmax = nmax, window = window, nanmath = nanmath)
    fconc && filter!(r -> r.Concentration != 0.0001, bold_df)
    logarithmic && (bold_df.Concentration = log10.(bold_df.Concentration))
    if method == :None
        gd = groupby(bold_df,[:Genotype, :Mouse, :Stim, :Area, :Concentration])
        pltdf = combine(gd, :Bold => (b -> mean(skipmissing(b))) => :Bold)
    else
        pltdf = transform(bold_df, :Bold => ByRow(b -> (ismissing(b) ? NaN : b)) => :Bold)
    end
    plt, plt_df, to_keep = prov_plot_bold(pltdf; nanmath = true, lim = lim)
    display(plot(plt..., size = (900*1.5,1200*1.5)))#layout = (6,3))
    ## find areas that respond to concentration
    m1, mp, to_keep2 = prov_odour_response(bold_df; farea = to_keep)
    # assess stim concentration interaction
    stim_df, m2, m3, p_complex, res = prov_stim_interaction(bold_df, to_keep2)
    show(m2)
    show(m3)
    show(p_complex)
    return bold_df, plt, plt_df, to_keep, m1, mp, to_keep2, stim_df, m2, m3, p_complex, res
end
##
bold_df, plt, plt_df, to_keep, m1, mp, to_keep2, stim_df, m2, m3, p_complex, res =
    prov_analysis_bold(fdata, :BOLD_90; fconc = true, av_pres = true, method = :Window,
    peak_to = :all, nmax = 6, window = 4:10, logarithmic = true,
    nanmath = true, lim = (-0.4,0.7))
##
union(bold_df.Genotype)
test_df = filter(r -> r.Genotype == "sert" && r.Area in to_keep && r.Area in to_keep2, bold_df)
f1 = @formula(Bold ~ 1 + Concentration * Area + Stim + Stim & Concentration + zerocorr(1+Concentration|Mouse) + zerocorr(1+Concentration|Area))
f2 = @formula(Bold ~ 1 + Concentration * Area * Stim + zerocorr(1+Concentration|Mouse) + zerocorr(1+Concentration|Area))
m1 = fit(MixedModel,f1,test_df)
m2 = fit(MixedModel,f2,test_df)
MixedModels.likelihoodratiotest(m1,m2)
##
f_add = @formula(Bold ~ 1 + Concentration * Area + Stim + Stim & Area + Stim & Concentration + zerocorr(1+Concentration|Mouse) + zerocorr(1+Concentration|Area))
f_mul = @formula(Bold ~ 1 + Concentration * Area + Stim + Stim & Concentration + Stim & Concentration & Area + zerocorr(1+Concentration|Mouse) + zerocorr(1+Concentration|Area))
m_add = fit(MixedModel,f_add,test_df)
m_mul = fit(MixedModel,f_mul,test_df)
##
function AICc(model)
    aic(model) + ((2*(dof(model)^2) + 2*dof(model))/(nobs(model) - dof(model) - 1))
end
function AICc_test(candidate, simpler)
    exp((AICc(candidate) - AICc(simpler))/2)
end

AICc_test(m_add,m_mul)


MixedModels.likelihoodratiotest(m_add,m2)
MixedModels.likelihoodratiotest(m_mul,m2)
##

df00 = fdata
signal = :BOLD_90
fconc = true
av_pres = true
method = :Window
peak_to = :all
nmax = 6
window = 4:10
logarithmic = true

df0 = prov_define_bold(df00, signal; av_pres = av_pres, method = :Window, window = window, nanmath = true)
fconc && filter!(r -> r.Concentration != 0.0001, df0)
logarithmic && (df0.Concentration = log10.(df0.Concentration))
open_html_table(df0)
pltdf = transform(df0, :Bold => ByRow(b -> (ismissing(b) ? NaN : b)) => :Bold)
plt, plt_df, to_keep = prov_plot_bold(pltdf)
plot(plt...)
sum(isnan.(df00.BOLD_90))
sum(ismissing.(df00.BOLD_90))
m1, mp, to_keep2 = prov_odour_response(df0; farea = to_keep)
to_keep2
