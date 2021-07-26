function xcur!()

    ax = zeros(mpnmax, mpnmax)

    mmaxp1 = Mmax + 1
    mmaxp2 = Mmax + 2

    bv = [-1., 0., 0.]

    for i in 1:Mmax
        for j in i:Mmax
            ax[i, j] = cl[i, j]
        end
        ax[i, Mmax + 1] = 0.
        for j in 1:Nmax
            ax[i, mmaxp1 + j] = bb[j, i]
        end
        if (icops >= 2)
            if(icops > 2)
                ax[i, mpnmax] = 0.
            else
                ax[i, mpnmax] = cl[mmaxp1, i]
            end
        end
    end

    ax[mmaxp1, mmaxp1] = 0.

    if (Nmax > 0)
        for j in 1:Nmax
            ax[mmaxp1, mmaxp1 + j] = bv[ityp[j]]
        end
    end
    if (icops >= 2)
        if (icops > 2)
            ax[mmaxp1, mpnmax] = 1.
        else
           ax[mmaxp1, mpnmax] = 0.
        end
    end

    for i in mmaxp2:mpnmax
        for j in i:mpnmax
            ax[i, j] = 0.
        end
    end

    for i in 1:mmaxp1
        fk[i] = 0.
    end

    if (llmax > 0)
        for ll in 1:llmax
            for i in 1:Mmax
                for j in i:Mmax
                    ax[i, j] = ax[i, j] + 2. * alph*eb[ll,i] * eb[ll, j]
                end
                ax[i, mmaxp1] = ax[i, mmaxp1] - 2. * alph*eb[ll,i]
                fk[i] += - 2. * alph*eb[ll,i] * eb[ll, mmaxp1]
            end
            ax[mmaxp1, mmaxp1] = ax[mmaxp1, mmaxp1] - 2. * alph
            fk[mmaxp1] += 2. * alph*eb[ll,mmaxp1]
        end
    end

    for i in 1:mpnmax
        for k in i:mpnmax
            ax[k, i] = ax[i, k]
        end
    end

    for j in 1:Nmax
        fk[mmaxp1 + j] = -bb[j, mmaxp1]
    end

    if (icops >= 2)
        fk[mpnmax] = value
    end

    gelg!(fk, ax, mpnmax, 1, 1e-7)

    global psicon = fk[mmaxp1]
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