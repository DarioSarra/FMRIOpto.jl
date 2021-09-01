include("prepareDB.jl")
##
fdata = filter(r ->
    r.Area != "PIR" &&
    r.Area != "OT",# &&
    #r.Valid_frs,
    # r.Area != "Ent" &&
    # r.Mouse != "FM14",
    fulldata)
    # lvA = union(fdata.Area)
    # fdata.Area = categorical(fdata.Area; levels = lvA)
## define bold
fconc = true
df0 = define_bold(fdata,:BOLD_85;av_pres = true, method = :Peak, peak_to = :all, nmax = 6, window = 4:12, nanmath = true)
# plot areas
    plt, meandf, to_keep = plot_bold(df0; nanmath = true, fconc = fconc)
    plot(plt..., size = (900*1.5,1200*1.5))#layout = (6,3)
## Areas that respond to concentration
m1, mp, to_keep2 = odour_response(df0; farea = to_keep, fconc = fconc)
    df1, m2, p_simple, p_complex, res = stim_interaction(df0, to_keep2; fconc = fconc)
        show(m2)
        show(p_simple)
        show(p_complex)
        res
##
bold_df, plt, plt_df, to_keep, m1, mp, to_keep2, stim_df, m2, p_simple, pcomplex, res =
    analysis_bold(fdata,:BOLD_90;
    method = :Window, av_pres = true, peak_to = :all,
    nmax = 6, window = 4:10, nanmath = true,
    fconc = true, lim = (-0.4,1))
    show(res)
#=
bold type
fconc = true use 3 highest concentration

method
    :Peak = from max onwards until peak_to & using the average response if av_pres = true
    :Max = which takes x nmax values from the average response per mouse and concentration
    :Window = average per animal of values in a window of frames (column Time)
