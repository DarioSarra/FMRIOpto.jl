function analyze_bold(df0, x)
    df = combine(groupby(df0,[:Area, :Stim, :Concentration]),  x .=> [mean, sem] .=>[:mean, :sem])
    df.Color = [x ? :red : :black for x in df.Stim]
    posdict = Dict(v=>p for (v,p) in zip(union(df.Concentration), collect(1:length(union(df.Concentration)))))
    df.Pos = [get(posdict, x,0) for x in df.Concentration]
    plt = []
        combine(groupby(df, :Area)) do dd
            pp = @df dd scatter(:Pos, :mean, yerror = :sem, group = :Stim, color = :Color, title = :Area[1], ylims = (-0.4,0.6), legend = false)
            @df dd plot!(:Pos, :mean, color = :Color, group = :Stim)
            Plots.abline!(0,0,color = :black, linestyle = :dash)
            push!(plt,pp)
        end
        pp = plot(plt..., layout = (6,3), size = (600*1.5,1200*1.5))

    df0.Test = df0[:,x]
    f = @formula(Test ~ 1 + Concentration*Area*Stim + (1|Mouse))
    # f = (MixedModels.term(x) ~ MixedModels.term(1) + MixedModels.term(:Concentration) &
    #     MixedModels.term(:Area) & MixedModels.term(:Stim) + (MixedModels.term(1) | MixedModels.term(:Mouse)))
    # terms = sum(MixedModels.term(t) for t in [1, :Concentration, :Area, :Stim])
    # inter = MixedModels.term(:Concentration) & MixedModels.term(:Area) & MixedModels.term(:Stim)
    # group = MixedModels.term(:Mouse)
    # response = MixedModels.term(x)
    # f = response ~ terms + (terms | group)
    m = fit(MixedModel, f,df0)
    f1 = @formula(Test ~ 1 + Concentration*Stim + (1|Mouse))
    mOFC = fit(MixedModel, f1,filter(r->r.Area == "OFC",df0))
    mOB = fit(MixedModel, f1,filter(r->r.Area == "OB",df0))
    maPIR = fit(MixedModel, f1,filter(r->r.Area == "aPIR",df0))
    return pp, df, m, mOFC, mOB, maPIR
end

function max_bold(v; max = 4)
    v2 = sort(v,rev = true)
    if length(v) <= max
        return mean(v2)
    else
        return mean(v2[1:max])
    end
end

function average_boldtime(b,t; window = 4:10)
    index = [x in window for x in t]
    mean(b[index])
end
