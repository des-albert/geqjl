function curf!()

    pp(x, alfa) = x^alfa * (2. - x^alfa)

    ff(x, alfa) = x^alfa * (2. - x^alfa)

    ppin(x, alfa) = (2. /(alfa + 1.)) * x^(alfa + 1.) - (1. /(2. * alfa + 1.)) * x^(2. * alfa + 1.)

    ffin(x, alfa) = (2. /(alfa + 1.)) * x^(alfa + 1.) - (1. /(2. * alfa + 1.)) * x^(2. * alfa + 1.)

    lcod = 0
    
    topol!(lcod)

    dint = zeros(10)
    beh = zeros(10)
    s = zeros(10)

    fmaxa = 0.
    mcont = 0
    rax = 0
    zax = 0

    for j in 1:Mr
        for n in 1:Nz
            if ( g[j,n] > fmaxa)
                fmaxa = g[j,n]
                jaxis = j
                naxis = n
                rax = R[jaxis]
                zax = Z[naxis]
            end
        end
    end

    global psicon
    global fabs = fmaxa + psicon

    if (mprfg)
        @printf(" Magnetic Axis radius = %12.5f height = %12.5f psi = %12.5f\n", rax,zax,fabs)
    end
    ipoi = 0
    idol = 0
    xalp = 0.

    for j = 1:Mr
        rzy = R[j]
        xalp1 = 0.
        for k in 1:10
            beh[k] = 0.
        end
        for n in 1: Nz
            xalp2 = 0.
            xalpr2 = 0.
            xalpz2 = 0.
            if ( g[j,n] > 0.)
                if ( ((j > 1) && (j < Mr)) && ((n > 1) && (n < Nz)) )
                    xalpr2 = ((g[j + 1, n] - g[j - 1, n]) / (2. * dr))^2
                    xalpz2 = ((g[j, n + 1] - g[j, n - 1]) / (2. * dz))^2

                    if ( g[j - 1, n] <= 0. )
                        xalpr2 = ((g[j + 1, n] - g[j,n])/dr)^2
                    end
                    if ( g[j + 1, n] <= 0. )
                        xalpr2 = ((g[j, n] - g[j - 1,n])/dr)^2
                    end
                    if ( g[j, n - 1] <= 0. )
                        xalpz2 = ((g[j, n + 1] - g[j,n])/dz)^2
                    end
                    if ( g[j, n + 1] <= 0. )
                        xalpz2 = ((g[j, n] - g[j,n - 1])/dz)^2
                    end

                    xalp2 += (xalpr2 + xalpz2)/rzy
                    ipoi += 1
                else
                    idol += 1
                end
                gn = g[j,n]/fmaxa
                s[1] = pp(gn, alfac) * rzy
                s[2] = ff(gn, alfac) / rzy
                s[3] = 1. /rzy^2
                s[4] = ppin(gn, alfac) * fmaxa
                s[5] = ffin(gn, alfac)/ rzy^2 * fmaxa
                s[6] = 1.
                s[7] = 2. * pi * rzy * ppin(gn, alfac)*fmaxa
                s[8] = 2. * pi * rzy
                s[9] = 0.
                s[10] = 0.

                for k in 1:10
                    beh[k] += s[k]
                end
            end
        end
        xalp += xalp1
        for k in 1:10
            dint[k] += beh[k]
        end


    end
    for k in 1:10
        dint[k] *= dz * dr
    end
    xalp = 2. * pi * xalp*dr*dz
    xind = 2. * xalp / (totcurr*totcurr*Rmpl)
    cjm = 0.
    c0pr = betapol*totcurr^2 / (8. * pi*dint[4])
    c0btr = (totcurr - c0pr*dint[1])/dint[2]

    if (mprfg)
        @printf(" c0pr = %12.5f c0btr = %12.5f\n",c0pr,c0btr)
        @printf(" Plasma Area = %12.5f  Plasma Volume = %12.5f\n\n",dint[6],dint[8])
    end

    for j in 1:Mr
        gn = g[j,naxis]/fmaxa
        pr[j] = c0pr*ppin(gn, alfac)*fmaxa
        bt2[j] = c0btr*ffin(gn, alfac)/R[j]^2 * fmaxa
    end

    for j in 1:Mr
        for n in 1:Nz
            gn = g[j,n]/fmaxa
            g[j,n] = (c0pr * pp(gn, alfac) * R[j]^2 + c0btr * ff(gn, alfac)) * dz^2
        end
    end

    pint = c0pr*dint[4]
    pintvo = c0pr*dint[7]
    bt2int = c0btr*dint[5]
    betap = 2. * pint/(totcurr^2 * dint[6])
    betapvo = 2. *pintvo/(totcurr^2 * dint[8])
    betat = dint[3]

    for j in 1:Mr
        cjt[j] = g[j,naxis]/R[j]
        if (abs(cjt[j]) > cjm)
            cjm = abs(cjt[j])
        end
    end
    if (cjm < dz^2)
        cjm = 1.
    end
    for j in 1:Mr
        cjt[j] = cjt[j]/cjm
    end

    return nothing

end
