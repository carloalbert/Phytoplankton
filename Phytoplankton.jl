function mle(data::AbstractVector{T}) where T <:AbstractFloat
    xmin = data[1]
    acc = zero(T)
    xlast = convert(T, Inf)
    ncount = 0
    for x in data
        if xlast == x
            continue
        end
        xlast = x
        ncount += 1
        acc += log(x / xmin)
    end
    ahat = 1 + ncount / acc
    stderr = (ahat - 1) / sqrt(ncount)
    return (ahat, stderr)
end
