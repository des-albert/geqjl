function bound_matrix!()
    ar = 2 / Mm1
    r0 = 1 / alpha
    az = sh*ar

    # Index vector interior neighbours to boundary, bottom-top

    j1 = Mr
    j2 = Mr * (Nz - 2)
    kp = 1

    for i in 2:Mm1
        ip[kp] = i + j1
        ip[kp + 1] = i + j2
        kp = kp + 2
    end

    # Left - Right

    for j in Mr:j2:Mr
        ip[kp] = j + 2
        ip[kp + 1] = j + Mm1
        kp += 2
    end
    kp = kp + 1

    # r - coordinates of boundary points, indexvector

    ra = r0 - 1
    rb = r0 + 1
    za = 0
    zb = Nm1*az
    m2 = Mr - 2
    j2 = Mr * Nm1

    rt = Vector{Float64}(undef, llp)
    zt = Vector{Float64}(undef, llp)

    # Bottom - Top

    nl = 1
    for i in 1:m2
        rt[nl] = ra + i*ar
        rt[nl + 1] = rt[nl]
        zt[nl] = za
        zt[nl + 1] = zb
        jp[nl] = i + 1
        jp[nl + 1] = 1 + 1 + j2
        nl += 2
    end
    
    # Left - Right

    n2 = Nz - 2
    for j in 1:n2
        rt[nl] = ra
        rt[nl + 1] = rb
        zt[nl] = j*az
        zt[nl + 1] = zt[nl]
        jp[nl] = j*Mr + 1
        jp[nl + 1] = (j + 1) * Mr
        nl += 2
    end

    # Matrix elements

    arh = ar/2.
    azh = az/2.

    for i in 1:kp

        nl = 1
        for j in 1:m2
            aux[nl, i] = flux(rt[i], rt[nl], zt[i] - zt[nl])/(sh*rt[nl])
            aux[nl + 1, i] = flux(rt[i], rt[nl + 1], zt[i] - zt[nl + 1])/(sh*rt[nl + 1])
            nl += 2
        end

        for j in 1:n2
            aux[nl, i] = flux(rt[i], rt[nl], zt[i] - zt[nl])*sh/(rt[nl] + arh)
            aux[nl + 1, i] = flux(rt[i], rt[nl + 1], zt[i] - zt[nl + 1])*sh/(rt[nl + 1] - arh)
            nl += 2
        end

    end

    return nothing
end

function flux(rv, rst, zv)

    ak = 4.0*rv*rst/((rv + rst)^2 + zv^2)
    x = ((rv - rst)^2 + zv^2)/((rv + rst)^2 + zv^2)
    return sqrt(rv*rst/ak)*((1 - ak/2.)*ellipk(x^2) - ellipe(x^2))/pi

end