function splnco!(psi)

    a = zeros(8)

    a1 = 4.
    for i in 1:4
        a[i] = 1. /a1
        a1 = a1 - 2.
        a[9-i] = 1. /a1
        a1 = a1*a1
    end

    # Control for inner loop on columns

    jh1 = 1
    lend = Mr*Nm1 + 1
    ls = Mr
    jh2 = min(16, Mm1 ÷ 2 )
    il = Mm1 - 1

@label  start
    k = 1
    jh = jh1
    mode = 2
@label middle
    j = 2 * jh

    for l in 1:ls:lend
        init = l + jh*mode
        iend = l + il
        for i in init:j:iend
            psi[i] = psi[i] + (psi[i + jh] + psi[i - jh] - 2. * psi[i]) * a[k]
        end
    end
    k += 1

    if (mode == 2)
        jh = 2 * jh
        if (jh < jh2)
            @goto middle
        end
        mode = 1
        if (k == 5)
            jh = jh ÷ 2
            if (jh >= jh1)
               @goto middle
            end
            if (jh1 == Mr)
               @goto last
            end
            jh1 = Mr
            lend = Mr
            ls = 1
            jh2 = min(Mr*16, Mr*(Nm1 ÷ 2))
            il = (Nm1 - 1)*Mr
            @goto start
        end
        k = 9 - k
        @goto middle

    elseif(mode == 1)
        jh = jh ÷ 2
        if(jh >= jh1)
            @goto middle
        end
        if(jh1 == Mr)
            @goto last
        end

        # Control for inner loop on rows

        jh1 = Mr
        lend = Mr
        ls = 1
        jh2 = min(Mr*16, Mr*(Nm1 ÷ 2))
        il = (Nm1 - 1)*Mr
        @goto start
    end
@label last
    return nothing

end