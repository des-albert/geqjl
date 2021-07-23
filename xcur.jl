function xcur!(a)

    mmaxp1 = Mmax + 1
    mmaxp2 = Mmax + 2

    bv = [-1., 0., 0.]

    for i in 1:Mmax
        for j in i:Mmax
            a[i, j] = cl[i, j]
        end
        a[i, Mmax + 1] = 0.
        for j in 1:Nmax
            a[i, mmaxp1 + j] = bb[j, i]
        end
        if (icops >= 2)
            if(icops > 2)
                a[i, mpnmax] = 0.
            else
                a[i, mpnmax] = cl[mmaxp1, i]
            end
        end
    end

    a[mmaxp1, mmaxp1] = 0.

    if (Nmax > 0)
        for j in 1:Nmax
            a[mmaxp1, mmaxp1 + j] = bv[ityp[j]]
        end
    end
    if (icops >= 2)
        if (icops > 2)
            a[mmaxp1, mpnmax] = 1.
        else
           a[mmaxp1, mpnmax] = 0.
        end
    end

    for i in mmaxp2:mpnmax
        for j in i:mpnmax
            a[i, j] = 0.
        end
    end

    for i in 1:mmaxp1
        fk[i] = 0.
    end

    if (llmax > 0)
        for ll in 1:llmax
            for i in 1:Mmax
                for j in i:Mmax
                    a[i, j] = a[i, j] + 2. * alph*eb[ll,i] * eb[ll, j]
                end
                a[i, mmaxp1] = a[i, mmaxp1] - 2. * alph*eb[ll,i]
                fk[i] = fk[i] - 2. * alph*eb[ll,i] * eb[ll, mmaxp1]
            end
            a[mmaxp1, mmaxp1] = a[mmaxp1, mmaxp1] - 2. * alph
            fk[mmaxp1] = fk[mmaxp1] + 2. * alph*eb[ll,mmaxp1]
        end
    end

    for i in 1:mpnmax
        for k in i:mpnmax
            a[k, i] = a[i, k]
        end
    end

    for j in 1:Nmax
        fk[mmaxp1 + j] = -bb[j, mmaxp1]
    end

    if (icops >= 2)
        fk[mpnmax] = value
    end

    gelg!(fk, a, mpnmax, 1, 1e-7)

    psicon = fk[mmaxp1]
    energy = 0.

    for i in 1:Mmax
        ki = i + 1
        if (ki <= Mmax)
            for k in ki:Mmax
                energy += cl[i, k] * fk[i] * fk[k]
            end
        end
    end

    fk[mmaxp2] = energy

    return nothing

end