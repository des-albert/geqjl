function splnco!(pss)

    local as = zeros(8)

    a1 = -4.
    for i in 1:4
        as[i] = 1. / a1
        a1 = a1 - 2.
        as[9 - i] = 1. / a1
        a1 = a1*a1
    end

    # Control for inner loop on columns

    jh1 = 1
    lend = Mr*Nm1 + 1
    ls = Mr
    jh2 = min(16, Mm1 ÷ 2 )
    il = Mm1 - 1

@label  L20
    k = 1
    jh = jh1
    mode = 2
@label L30
    j = 2 * jh

    for l in 1:ls:lend
        init = l + jh*mode
        iend = l + il
        for i in init:j:iend
            pss[i] = pss[i] + (pss[i + jh] + pss[i - jh] - 2. * pss[i]) * as[k]
        end
    end
    k += 1

    if (mode == 2)
        jh = 2 * jh
        if (jh < jh2)
            @goto L30
        end
        mode = 1
        if (k == 5)
            jh = jh ÷ 2
            if (jh >= jh1)
               @goto L30
            end
            if (jh1 == Mr)
               @goto L50
            end
            jh1 = Mr
            lend = Mr
            ls = 1
            jh2 = min(Mr * 16, Mr * (Nm1 ÷ 2))
            il = (Nm1 - 1) * Mr
            @goto L20
        end
        k = 9 - k
        @goto L30

    elseif (mode == 1)

        jh = jh ÷ 2
        if(jh >= jh1)
            @goto L30
        end
        if(jh1 == Mr)
            @goto L50
        end

        # Control for inner loop on rows

        jh1 = Mr
        lend = Mr
        ls = 1
        jh2 = min(Mr * 16, Mr * (Nm1 ÷ 2))
        il = (Nm1 - 1) * Mr
        @goto L20
    end
@label L50
    return nothing

end