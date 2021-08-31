nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,y)
nansem(x) = sem(filter(!isnan,x))
nansem(x,y) = mapslices(nansem,x,y)

function safemean(x)
    if sum(.!(isnan.(x))) < 4
        return NaN
    else
        return nanmean(x)
    end
end

function safesem(x)
    if sum(.!(isnan.(x))) < 4
        return NaN
    else
        return nansem(x)
    end
end
##
"""
   definBold(df, signal; method = :Max, pick_to = :all, nmax = 4, window = 4:10, nanmath = true)
It adjust the Bold in column signal per genotype, mouse, stim, area and concentration
according to different algorythms define by method options listed below
    :Pick = take all values after the pick until the index define by pick_to
        if pick_to == :all takes all values available
    :Max = take the highest n values defined by nmax, default = 4
    :Window = averages over the time points specified by window, default = 4:10
nanmath: if true (default) uses safemean and safe sem to calculate mean and sem,
else it uses regular mean and sem

"""

function defineBold(df, signal; method = :Pick, pick_to = :all, nmax = 4, window = 4:10, nanmath = true)
    fmean, fsem = nanmath ? (safemean, safesem) : (mean, sem)
    gd0 = groupby(df,[:Genotype, :Mouse, :Stim, :Area, :Concentration])
    if method == :Pick
        df0 = pick_bold(gd0, signal; pick_to = pick_to)
        df0.Bold = df0[:,signal]
    elseif method == :Max
        df0 = combine(gd0, signal => (s -> max_bold(s; max = nmax, nanmath = nanmath)) => :Bold)
    elseif method == :Window
        df0 = combine(gd0, [signal, :Time] =>((s,t) -> window_bold(s,t; window = window, nanmath = nanmath)) => :Bold)
    else
        error("unrecognised bold processing method:/n options are :Pick, :Max, :Timewindow")
    end
    return df0
end

function pick_bold(gd, signal; pick_to = :all)
    combine(gd) do dd
        any(ismissing,dd[:,signal]) && println("got NaN ", dd.Mouse[1], " ", dd.Area[1], " ", dd.Stim[1], " ", dd.Concentration[1])
        if all(ismissing.(dd[:,signal]))
            return filter(r -> r.Area == "NoData",dd)
        else
            i = findmax(skipmissing(dd[:,signal]))[2]
            return pick_to == :all ? dd[i:end,:] : dd[i:pick_to,:]
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

function plot_bold(df0,x; nanmath = true)
    df1 = combine(groupby(df0,[:Mouse,:Area, :Stim, :Concentration]),  x .=> (t -> mean(skipmissing(t))) .=> x)
    if nanmath
        df = combine(groupby(df1,[:Area, :Stim, :Concentration]),  x .=> [safemean, safesem] .=>[:mean, :sem])
    else
        df = combine(groupby(df1,[:Area, :Stim, :Concentration]),  x .=> [mean, sem] .=>[:mean, :sem])
    end
    df.Color = [x ? :red : :black for x in df.Stim]
    posdict = Dict(v=>p for (v,p) in zip(union(df.Concentration), collect(1:length(union(df.Concentration)))))
    df.Pos = [get(posdict, x,0) for x in df.Concentration]
    plt = []
        combine(groupby(df, :Area)) do dd
            pp = @df dd scatter(:Pos, :mean, yerror = :sem, group = :Stim, color = :Color,
                title = :Area[1], legend = false)
            @df dd plot!(:Pos, :mean, color = :Color, group = :Stim)
            Plots.abline!(0,0,color = :black, linestyle = :dash)
            push!(plt,pp)
        end
        #pp = plot(plt..., size = (600*1.5,1200*1.5))#layout = (6,3)
    return plt, df
end
