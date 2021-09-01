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

iscolumn(df,x) = x in propertynames(df)

function model_res(m)
    DataFrame(Area = coefnames(m), P = m.pvalues, Z = m.beta ./ m.stderror)
end
