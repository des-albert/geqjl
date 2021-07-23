function gelg!(rr, a, m, n, eps)

    if (m > 0)
        ier = 0
        piv = 0.
        mm = m*m
        nm = n*m
        for l in 1:mm
            tb = abs(a[l])
            if (tb > piv)
                piv = tb
                i = l
            end
        end
        tol = eps*piv

        # a(i) is pivot element. piv contains the absolute value of a(i)

        # Start elimination loop

        lst = 1
        for k in 1:m

        # Test on singularity

            if (piv <= 0)
                @goto last
            end
            if (ier == 0)
                if (piv <= tol)
                    ier = k - 1
                end
            end
            pivi = 1. / a[i]
            j = (i - 1) ÷ m
            i = i - j*m - k
            j = j + 1 - k

            # i+k is row-index, j+k column-index of pivot element
            # pivot row reduction and row interchange in right hand side r

            for l in k:m:nm
                ll = l + i
                tb = pivi*rr[ll]
                rr[ll] = rr[l]
                rr[l] = tb
            end

            # Is elimination terminated

            if (k >= m)
                @goto eliminend
            end

            # Column interchange in matrix a

            lend = lst + m - k
            if (j > 0)
                ii = j*m
                for l in lst:lend
                    tb = a[l]
                    ll = l + ii
                    a[l] = a[ll]
                    a[ll] = tb
                end
            end

            #  Row interchange and pivot row reduction in matrix a

            for l in lst:m:mm
                ll = l + i
                tb = pivi*a[ll]
                a[ll] = a[l]
                a[l] = tb
            end

            # Save column interchange information

            a[lst] = j

            # Element reduction and next pivot search
!
            piv = 0.
            lst += 1
            j = 0
            for ii in lst:lend
                pivi = -a[ii]
                ist = ii + m
                j += 1
                for l in ist:m:mm
                    ll = l - j
                    a[l] += pivi*a[ll]
                    tb = abs(a[l])
                    if (tb > piv)
                        piv = tb
                        i = l
                    end
                end
                for l in k:m:nm
                    ll = l + j
                    rr[ll] = rr[ll] + pivi*rr[l]
                end
            end
            lst += m

        end

@label eliminend

        # End of elimination loop
        # Back substitution and back interchange

        if (m >= 1)
            if (m != 1)
                ist = mm + m
                lst = m + 1
                for i in 2:m
                    ii = lst - i
                    ist = ist - lst
                    l = ist - m
                    l = floor(Int, a[l] + 0.5)
                    for j in ii:m:nm
                        tb = rr[j]
                        ll = j
                        for k in ist:m:mm
                            ll += 1
                            tb = tb - a[k] * rr[ll]
                        end
                        k = j + l
                        rr[j] = rr[k]
                        rr[k] = tb
                    end
                end
            end
            return nothing
        end

    end
@label  last
    ier = -1
    return nothing
end