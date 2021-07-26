function topol!(lcod)

    ndi = [1, -1]

    lsym = 1
    nmx = zeros(2)
    if (naxis != 1)
        lsym = 2
    end
    lcod = 0
    js = 0
    global psicon

    for j in 1:Mr
        for n in 1:Nz
            g[j, n] -= psicon
            if (g[j, n] <= 0.)
                g[j, n] = 0.
            end
        end
    end

    n = naxis
    jma = jaxis
    if (g[jma, n] == 0.)
        lcod = 5
        println(" Plasma disappeared of region, probably squeezed of f at center ")
    else
        for l in 1:lsym
            jma = jaxis
            n = naxis
            for j in jma:Mr
                if (g[j, n] == 0.)
                    js = j
                    @goto L10
                end
            end
            @goto L100
@label L10
            jmax = js - 1
            jnu = Mr - jmax + 1
            for ji in jnu:Mr
                j = Mr - ji + 1
                if (g[j, n] == 0.)
                    js = j
                    @goto L20
                end
            end
            @goto L110

@label L20
            jmin = js + 1
            while(true)
                j1 = jmin - 1
                j2 = jmax + 1
                for j in 1:j1
                    g[j, n] = 0.
                end
                for j in j2:Mr
                    g[j, n] = 0.
                end
                if ( jmax <= jmin)
                    @goto L90
                end
                n += ndi[l]
                if ( ! ((n < Nz) && (n > 1)))
                    @goto L140
                end
                if (g[jmax,n] == 0.)
                    jnu = Mr - jmax + 1
                    for ji in jnu:Mr
                        j = Mr - ji + 1
                        if ( g[j,n] != 0.)
                            js = j
                            @goto L30
                        end
                    end
                    @goto L70
@label L30
                    jmax = js
                else
                    jma = jmax
                    for j in jma:Mr
                        if ( g[j,n] == 0.)
                            js = j
                            @goto L40
                        end
                    end
                    @goto L120

@label L40
                    jmax = js - 1
                end
                jmi = jmin
                if (g[jmin, n] == 0.)
                    for j in jmi:Mr
                        if ( g[j,n] != 0.)
                            js = j
                            @goto L50
                        end
                    end
                    @goto L80
@label L50
                    jmin = js
                else
                    jnu = Mr - jmin + 1
                    for ji in jnu:Mr
                        j = Mr - ji + 1
                        if (g[j, n] == 0.)
                            js = j
                            @goto L60
                        end
                    end
                    @goto L130
@label L60
                    jmin = js + 1
                end
            end
@label L70
            n -= ndi[l]
            @goto L90
@label L80
            n -= 1
@label L90
            nmx[l] = n + ndi[l]
        end

        for n in 1:Nz
            if ( (nmx[1] - n)*(nmx[2] - n) >= 0. )
                for j= 1:Mr
                    g[j,n] = 0.
                end
            end
        end

        return nothing

@label L100
        println(" Plasma runs out of grid on outside at axis")
        lcod = 1
        @goto L150

@label L110
        println(" Plasma runs out of grid on inside at axis")
        lcod = 6
        @goto L150

@label L120
        println("'Plasma runs out of grid on outside for n = '", n)
        lcod = 3
        @goto L150

@label L130
        println("'Plasma runs out of grid on inside for n = '", n)
        lcod = 4
        @goto L150

@label L140
        println("Plasma runs out on top immersing probably top conductor")
        lcod = 2
        @goto L150
    end
@label L150
    println(" Case abondoned")
    return nothing
end