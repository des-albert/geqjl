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

    for j in Mr:Mr:j2
        ip[kp] = j + 2
        ip[kp + 1] = j + Mm1
        kp += 2
    end
    kp -= 1

    # r - coordinates of boundary points, indexvector

    ra = r0 - 1.
    rb = r0 + 1.
    za = 0.
    zb = Nm1*az
    m2 = Mr - 2
    j2 = Mr * Nm1

    rt = zeros(llp)
    zt = zeros(llp)

    # Bottom - Top

    nl = 1
    for i in 1:m2
        rt[nl] = ra + i*ar
        rt[nl + 1] = rt[nl]
        zt[nl] = za
        zt[nl + 1] = zb
        jp[nl] = i + 1
        jp[nl + 1] = i + 1 + j2
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
            aux[nl, i] = gfl(rt[i], rt[nl], zt[i] - zt[nl], ar)/(sh*rt[nl])
            aux[nl + 1, i] = gfl(rt[i], rt[nl + 1], zt[i] - zt[nl + 1], ar)/(sh*rt[nl + 1])
            nl += 2
        end

        for j in 1:n2
            aux[nl, i] = gfl(rt[i], rt[nl], zt[i] - zt[nl], az)*sh/(rt[nl] + arh)
            aux[nl + 1, i] = gfl(rt[i], rt[nl + 1], zt[i] - zt[nl + 1], az)*sh/(rt[nl + 1] - arh)
            nl += 2
        end

    end

    return nothing
end

function gfl(rv, rst, zv, del)    
    ak = 4.0*rv*rst/((rv + rst)^2 + zv^2)   
    y = ((rv - rst)^2 + zv^2)/((rv + rst)^2 + zv^2)
    if (y == 0.)
        xdl = 2. * (log(del/(4. *rv)) - 1.)
    else
        xdl = log(y)
    end

    return sqrt(rv*rst/ak)*((1. - ak/2.)*elk(y, xdl) - ele(y, xdl))/pi
end

    p1(x) = (((.01736506451*x + .04757383546)*x + .06260601220)*x + .44325141463)*x + 1.0
    p2(x) = (((.00526449639*x + .04069697526)*x + .09200180037)*x + .24998368310)*x
    p3(x) = (((.01451196212*x + .03742563713)*x + .03590092383)*x + .09666344259)*x + 1.38629436112
    p4(x) = (((.00441787012*x + .03328355346)*x + .06880248576)*x + .12498593597)*x + .5
    elk(x, xdl) = p3(x) - xdl*p4(x)
    ele(x, xdl) = p1(x) - xdl*p2(x)

function fl(rv, rst, zv)

    ak = 4.0*rv*rst/((rv + rst)^2 + zv^2)
    x = clamp(((rv - rst)^2 + zv^2)/((rv + rst)^2 + zv^2), 5.e-6, 1.0)
    return sqrt(rv*rst/ak)*((1. - ak/2.)*ellipk(1. - x) - ellipe(1. - x))/pi

end