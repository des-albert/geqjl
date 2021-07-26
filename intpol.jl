function condit!(vv, rx, zx, itp, xv, i, j)
    xc = (rx - Rmin) / dr
    yc = (zx - Zmin) / dz

    var = zeros(3)
    intpol!(vv, xc, yc, var)
    xv[i, j] = var[itp]

    return nothing
end

function intpol!(tt, xi, yi, w)

    p(c) = (a1*c + a2) * c / 2. + 0.125 * a1

    dp(c) = a1*c + a2 / 2.

    af = zeros(3)
    adf = zeros(3)

    dx = xi - floor(Int, xi + 0.5)
    dy = yi - floor(Int, yi + 0.5)

    m = floor(Int, xi - 0.5)
    n = floor(Int, yi + 0.5)

    mni = m + n*Mr

    if (! (m < 0 || m >= (Mr - 2) || n < 0 || n >= (Nz - 1)))
        if (n == 0)
            for i in 1:3
                ki = mni + 1
                a1 = 2. * (tt[k + Mr] - tt[ki])
                a2 = 0.
                af[i] = p(dy) + tt[ki]
                adf[i] = dp(dy)
           end
        else
            for i in 1:3
                ki = mni + i
                a1 = tt[ki + Mr] + tt[ki - Mr] - 2. * tt[ki]
                a2 = tt[ki + Mr] - tt[ki - Mr]
                af[i] = p(dy) + tt[ki]
                adf[i] = dp(dy)
            end
        end

        a1 = af[1] + af[3] - 2. * af[2]
        a2 = af[3] - af[1]
        w[1] = p(dx) + af[2]
        w[3] = dp(dx)
        a1 = adf[1] + adf[3] - 2. * adf[2]
        a2 = adf[3] - adf[1]
        w[2] = p(dx) + adf[2]

    end
    return nothing
end

