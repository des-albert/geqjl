function condit!(v, rx, zx, itp, b)
    x = (rx - Rmin)/dr
    y = (zx - Zmin)/dz

    var = zeros(3)
    intpol!(v, x, y, var[1], var[3], var[2])
    b = var[itp]

    return nothing
end

function intpol!(v, x, y, fi, fx, fy)

    p(x) = (a1*x + a2) * x/2. + 0.125 * a1

    dp(x) = a1*x + a2/2.

    af = zeros(3)
    adf = zeros(3)

    dx = x - floor(x + 0.5)
    dy = y - floor(y + 0.5)

    m = floor(Int, x - 0.5)
    n = floor(Int, y + 0.5)

    mni = m + n*Mr

    if (! (m < 0 || m >= (Mr - 2) || n < 0 || n >= (Nz - 1)))
        if (n == 0)
            for i in 1:3
                ki = mni + 1
                a1 = 2. * (v[k + Mr] - v[ki])
                a2 = 0.
                af[i] = p(dy) + v[ki]
                adf[i] = dp(dy)
           end
        else
            for i in 1:3
                ki = mni + i
                a1 = v[i + Mr] + v[ki - Mr] - 2. * v[ki]
                a2 = v[ki + Mr] - v[ki - Mr]
                af[i] = p(dy) + v[ki]
                adf[i] = dp(dy)
            end
        end

        a1 = af[1] + af[3] - 2. * af[2]
        a2 = af[3] - af[1]
        fi = p(dx) + af[2]
        fx = dp(dx)
        a1 = adf[1] + adf[3] - 2. * adf[2]
        a2 = adf[3] - adf[1]
        fy = p(dx) + adf[2]

    end
    return nothing
end

