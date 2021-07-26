function flux!(psi)
 
    local tc = zeros(llp)
    local aa = zeros(llp)

    nn = Nm1
    iu = Mm1 - 1
    itc = 2 * Mm1
    lo = nn ÷ 2
    tc[itc + lo] = 0.

@label L10
    l = lo ÷ 2
    tc[itc + l] = sqrt(2. + tc[itc + lo])
    lo = l
@label L20
    tc[itc + nn - l] = -tc[itc + l]
    l = l + 2 * lo

    sw = ((2 * l ÷ nn) * (2 * lo - 3))

    if (sw > 0)
        @goto L10
    elseif (sw == 0)
        tc[itc + l] = (tc[itc + l + lo] + tc[itc + l - lo]) / tc[itc + lo]
        @goto L20
    elseif (sw < 0)
       for i in 1:iu
            d = alpha / (Mm1 + alpha * (2 * i - Mm1))
            tc[i] = ss / (1. - d)
            tc[i + Mm1] = ss / (1. + d)
        end
    end

    lo = nn ÷ 2
    j1 = 1 + Mr * (Nm1 ÷ nn)
    ju = Nm1*Mr

    for i in j1:Mr:ju
        psi[i + 1] = psi[i + 1] + tc[1] * psi[i]
        psi[iu + i] = psi[iu + i] + tc[iu + Mm1] * psi[i + Mm1]
    end

    aa[Mm1] = 0.
    aa[Mm1 + Mm1] = 0.
    mode = 2
    is = -1
@label  L80
    li = 2 * lo

    iphase = 2 * mode - (li ÷ nn)
    jd = Mr * nn ÷ li
    jh = jd ÷ 2
    jt = jd + jh
    ji = 2 * jd
    jo = jd * mode * ((1 - is) ÷ 2) + 1

    for j in jo:ji:ju
        j1 = j + 1
        jdm = jd * is
        jhm = jh * is
        jtm = jt * is
        jiu = j + iu

        if (iphase == 1)
            for i in j1:jiu
                aa[i - j] = psi[i] + psi[i + jd] + psi[i + jdm]
                psi[i] = 0.0
            end
        elseif (iphase == 2)
            for i in j1:jiu
                aa[i - j] = psi[i] + psi[i + jd] + psi[i + jdm]
                psi[i] = (psi[i] - psi[i + jh] - psi[i + jhm]) / 2.
            end
        elseif (iphase == 3)
            for i in j1:jiu
                aa[i - j] = 2. * psi[i]
                psi[i] = psi[i + jd] + psi[i + jdm]
            end
        elseif (iphase == 4)
            for i in j1:jiu
                d = psi[i] - psi[i + jt] - psi[i + jtm]
                psi[i] = psi[i] - psi[i + jh] - psi[i + jhm] + psi[i + jd] + psi[i + jdm]
                aa[i - j] = d + psi[i]
            end
        end

        for l in lo:li:nn
            d = 2. - tc[itc + l]
            for i in 1:iu
                k = Mm1 - i
                b = 1. / (d + tc[k] + tc[k + Mm1] * (1. - aa[k + Mr]))
                aa[k + Mm1] = tc[k] * b
                aa[k] = (aa[k] + tc[k + Mm1] * aa[k + 1]) * b
            end
            for i in 2:iu
                aa[i] = aa[i + Mm1] * aa[i - 1] + aa[i]
            end
        end
        for i in j1:jiu
            psi[i] += aa[i - j]
        end
        is = -1

    end


    if (iphase == 1)
        return nothing
    elseif (iphase == 2)
        lo = 2 * lo
        if (lo < nn)
            @goto L80
        end
    elseif (iphase == 3 || iphase == 4)
        lo = lo ÷ 2
        if (lo == 1)
            mode = 1    
        end
        @goto L80
    end

    return nothing
end