"""
   defineBold(df, signal; method = :Max, peak_to = :all, nmax = 4, window = 4:10, nanmath = true)
It adjust the Bold signal per genotype, mouse, stim, area and concentration
according to different algorythms define by method option.
:Peak = take all values after the peak until the index define by peak_to
    if peak_to == :all takes all values available (goes over every presentation
    if the presentation where not averaged before)
:Max = take the highest n values defined by nmax, default = 4
:Window = averages over the time points specified by window, default = 4:10
nanmath: if true (default) uses safemean and safe sem to calculate mean and sem,
else it uses regular mean and sem

"""
function define_bold(df0, signal; av_pres = false, method = :Peak, peak_to = :all, nmax = 4, window = 4:10, nanmath = true)
    fmean, fsem = nanmath ? (safemean, safesem) : (mean, sem)
    if method == :Peak
        if av_pres
            #average response over presentation
            df1 = average_presentation(df0,signal; nanmath = nanmath)
            signal = :Bold
            gd1 = groupby(df1,[:Genotype, :Mouse, :Stim, :Area, :Concentration])
        else
            #change groupby to find peak for each presentation
            gd1 = groupby(df0,[:Genotype, :Mouse, :Stim, :Area, :Concentration,:Repetition_nr])
        end
        #find peak in the bold signal and take all value from there until peak_to
        df2 = peak_bold(gd1, signal; peak_to = peak_to)
        df2.Bold = df2[:,signal]
    elseif method == :Max
        #average response over all presentation
        df1 = average_presentation(df0,signal; nanmath = nanmath)
        gd1 = groupby(df1,[:Genotype, :Mouse, :Stim, :Area, :Concentration])
        signal = :Bold
        #find the first n max values in the bold signal
        df2 = combine(gd1, signal => (s -> max_bold(s; max = nmax, nanmath = nanmath)) => :Bold)
    elseif method == :Window
        #average all bold values per mouse and concentration in a time window
        df1 = average_presentation(df0,signal; nanmath = nanmath)
        gd1 = groupby(df1,[:Genotype, :Mouse, :Stim, :Area, :Concentration])
        signal = :Bold
        df2 = combine(gd1, [signal, :Time] =>((s,t) -> window_bold(s,t; window = window, nanmath = nanmath)) => :Bold)
    else
        error("unrecognised bold processing method:/n options are :Peak, :Max, :Timewindow")
    end
    return df2
end

function average_presentation(df0,signal; nanmath = true)
    fmean, fsem = nanmath ? (safemean, safesem) : (mean, sem)
    df0[:,signal] = [ismissing(x) ? NaN : x for x in df0[:,signal]]
    df1 = combine(groupby(df0, [:Time,:Concentration,:Area,:Stim,:Mouse,:Genotype]), signal => fmean => :Bold)
    return df1
end

function peak_bold(gd, signal; peak_to = :all)
    combine(gd) do dd
        #any(ismissing,dd[:,signal]) && println("got NaN ", dd.Mouse[1], " ", dd.Area[1], " ", dd.Stim[1], " ", dd.Concentration[1])
        if all(ismissing.(dd[:,signal]))
            return filter(r -> r.Area == "NoData",dd)
        else
            i = findmax(skipmissing(dd[:,signal]))[2]
            return peak_to == :all ? dd[i:end,:] : dd[i:peak_to,:]
        end
    end
end

function max_bold(v; max = 4, nanmath = true)
    fmean = nanmath ? safemean : mean
    v1 = [ismissing(x) ? NaN : x for x in v]
    v2 = sort(v1,rev = true)
    if length(v1) <= max
        return fmean(v2)
    else
        return fmean(v2[1:max])
    end
end

function window_bold(b,t; window = 4:10, nanmath = true)
    fmean = nanmath ? safemean : mean
    index = [x in window for x in t]
    b1 = [ismissing(x) ? NaN : x for x in b]
    fmean(b1[index])
end

"""
    plot_bold(df0,x; nanmath = true)
Function to check the average response over concentration per area.
df0 has to be a dataframe preprocessed with define_bold,
containing bold signal in a column named :Bold.
returns a vector of plots of bold over concentration per each area,
a dataframe with the mean and sem values of bold per area
and a filtered list of the areas that contain sufficient data for further analysis
"""

function plot_bold(df0; nanmath = true, fconc = true, lim = (-0.4,0.7))
    fmean = nanmath ? safemean : mean
    fsem = nanmath ? safesem : sem
    #calculate mean repsonse over concentrations per area
    df = combine(groupby(df0,[:Area, :Stim, :Concentration, :Genotype]),  :Bold .=> [fmean, fsem] .=>[:mean, :sem])
    # identify regions with missing concentration info
    res = combine(groupby(df,[:Area]), :mean => (m-> .!(any(isnan.(m)))) => :to_keep)
    to_keep = filter(r -> r.to_keep, res)[:,:Area]
    # prepare a grid plot with the results for each regions
    df.Color = [x ? :red : :black for x in df.Stim]
    posdict = Dict(v=>p for (v,p) in zip(union(df.Concentration), collect(1:length(union(df.Concentration)))))
    df.Pos = [get(posdict, x,0) for x in df.Concentration]
    sort!(df,[:Stim,:Pos])
    fconc && filter!(r -> r.Concentration != 0.0001 ,df)
    df1 = filter(r-> r.Genotype == "sert",df)
    plt = []
        combine(groupby(df1, :Area)) do dd
            pp = @df dd scatter(:Pos, :mean, yerror = :sem, group = :Stim, color = :Color,
                title = :Area[1], legend = false, ylims = lim)
            @df dd plot!(:Pos, :mean, color = :Color, group = :Stim)
            Plots.abline!(0,0,color = :black, linestyle = :dash)
            push!(plt,pp)
        end
    return plt, df, to_keep
end

"""
    odour_response(df0, to_keep)
Use control data to identify which areas show a dose dependent response to
odour concentration with a mixed model.
df0 has to be a dataframe preprocessed with define_bold,
    containing bold signal in a column named :Bold.
filt is a vector of string specifying which area to keep and which to filter
    out from the analysis, default is :none which takes all areas
"""

function odour_response(df0; farea = :none, fconc = true)
    # use only control data to determine odour concentration dose dependency
    df1 = filter(r -> !r.Stim, df0)
    #remove lowest concentration
    fconc && filter!(r -> r.Concentration != 0.0001, df1)
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
    #replace NaN values used for plotting to missing for model
    df1.Bold = [isnan(x) ? missing : x for x in df1.Bold]
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

"""
    stim_interaction(df0, farea; fconc = :none)
Use  data to identify which areas show an interaction between dose dependent
response to odour concentration and stimulation with a mixed model.
df0 has to be a dataframe preprocessed with define_bold,
    containing bold signal in a column named :Bold.
farea is a vector of string specifying which areas showed a significant dose
dependent response to odour concetration, found using odour_response function
fconc allows to remove lowest concetration from the dataset
"""
function stim_interaction(df0, farea; fconc = true)
    # filters area lacking concentration response identified by odour_response function
    df1 = filter(r -> r.Genotype == "sert" && r.Area in farea, df0)
    fconc && filter!(r -> r.Concentration != 0.0001, df0)
    # identify if the model requires to incorporate time
    if iscolumn(df0,:Time)
        f1 = @formula(Bold ~ 1 + Time + Concentration*Area + Stim + (1+Time|Mouse))
        f2 = @formula(Bold ~ 1 + Time + Concentration + Area + Stim + Concentration & Area + Concentration & Area & Stim + (1+Time|Mouse))
        f3 = @formula(Bold ~ 1 + Time + Concentration*Area*Stim + (1+Time|Mouse))
    else
        f1 = @formula(Bold ~ 1 + Concentration*Area + Stim + (1|Mouse))
        f2 = @formula(Bold ~ 1 + Concentration + Area + Stim + Concentration & Area + Concentration & Area & Stim + (1|Mouse))
        f3 = @formula(Bold ~ 1 + Concentration*Area*Stim + (1|Mouse))
    end
    #replace NaN values used for plotting to missing for model
    df1.Bold = [isnan(x) ? missing : x for x in df1.Bold]
    m1 = fit(MixedModel, f1,df1)
    m2 = fit(MixedModel, f2,df1)
    m3 = fit(MixedModel, f3,df1)
    p_simple = MixedModels.likelihoodratiotest(m1,m2)
    p_complex = MixedModels.likelihoodratiotest(m2,m3)
    res = model_res(m2)
    filter!(r -> occursin("Stim", r.Area), res)
    sort!(res,:P)
    return df1, m2, p_simple, p_complex, res
end

function analysis_bold(df0, signal; fconc = true, av_pres = true, method = :Peak, peak_to = :all, nmax = 6, window = 4:12, nanmath = true, lim = (-0.4,0.7))
    #adjust bold signal
    bold_df = define_bold(df0,signal; method = method, av_pres = av_pres, peak_to = peak_to, nmax = nmax, window = window, nanmath = nanmath)
    plt, plt_df, to_keep = plot_bold(bold_df; nanmath = true, fconc = fconc, lim = lim)
    display(plot(plt..., size = (900*1.5,1200*1.5)))#layout = (6,3))
    ## find areas that respond to concentration
    m1, mp, to_keep2 = odour_response(bold_df; farea = to_keep, fconc = fconc)
    # assess stim concentration interaction
    stim_df, m2, p_simple, p_complex, res = stim_interaction(bold_df, to_keep2; fconc = fconc)
    show(m2)
    show(p_simple)
    show(p_complex)
    return bold_df, plt, plt_df, to_keep, m1, mp, to_keep2, stim_df, m2, p_simple, p_complex, res
end
