module geq
    using Base:Float64  
    using DelimitedFiles
    using Printf

    include("bound_matrix.jl")
    include("eqsil.jl")
    include("flux.jl")
    include("splnco.jl")
    include("intpol.jl")
    include("startt.jl")
    include("compar.jl")
    include("topol.jl")
    include("curf.jl")
    include("saddle.jl")
    include("xcur.jl")
    include("gelg.jl")
    include("plotit.jl")

    Nexp = 6
    Mc = 15
    Ng = 1
    Nmax = 6
    llmax = 16

    Mr = 1 + 2^Nexp
    Nz = Mr
    MN = Mr * Nz
    Mm1 = Mr - 1
    Nm1 = Nz - 1
    llp = 2 * (Mr + Nz) - 8


    dataFile = open("solver.dat", "r")
    Rmin, Rmax, Zmin, Zmax, error, meshfg, mprfg = readdlm(IOBuffer(readline(dataFile)))
    Rmpl, Offset, Apl, El, tri, Rxpsn, Zxpsn = readdlm(IOBuffer(readline(dataFile)))

    if (meshfg)
        Rmin = Rmpl - 32 * Apl / 20
        Rmax = Rmpl + 32 * Apl / 20
        Zmin = Offset - 32 * (Offset + abs(Zxpsn)) / 20
        Zmax = Offset + 32 * (Apl * El) / 20
    end


    ip = zeros(Int, llp)
    jp = zeros(Int, llp)
    aux = zeros(llp, llp)

    dr = (Rmax - Rmin) / Mm1
    dz = (Zmax - Zmin) / Nm1
    alpha = (Rmax - Rmin) / (Rmax + Rmin)
    sh = dz / dr
    ss = sh^2

    R = zeros(Mr)
    Z = zeros(Nz)
    cjt = zeros(Mr)
    pr = zeros(Mr)
    bt2 = zeros(Mr)

    for i in 1:Mr
        R[i] = Rmin + (i - 1) * dr
    end

    for j in 1:Nz
        Z[j] = Zmin + (j - 1) * dz
    end

    @printf("Rmin  = %7.3f Rmax = %7.3f dr = %7.3f\n", Rmin, Rmax, dr)
    @printf("Zmin  = %7.3f Zmax = %7.3f dz = %7.3f\n", Zmin, Zmax, dz)
    @printf("alpha = %7.3f sh   = %7.3f ss = %7.3f\n", alpha, sh, ss)

    bound_matrix!()

    # Read Poloidal Field Coil Data

    ic = zeros(Int, 15)
    Ra = zeros(Ng, Mc)
    Za = zeros(Ng, Mc)
    Ex = zeros(Ng, Mc)
    Rl = zeros(Ng, Mc)

    k = 0
    while (true)
        global k = k + 1
        ic[k] = 0
        while (true)
            rac, zac, exc = readdlm(IOBuffer(readline(dataFile)))
            if (rac > 0.0)
                ic[k] += 1
                Ra[ic[k], k] = rac
                Za[ic[k], k] = zac
                Ex[ic[k], k] = exc
                Rl[ic[k], k] = 1.e-20 * rac
            else
                break
            end
        end
        if (ic[k] < 1)
            break
        end
    end

    Mmax = k - 1

    if (mprfg)
        println(" Conductor groups available for optimization")
        for k in 1:Mmax
            println("Group : ", k)
            for i in 1:ic[k]
                @printf("  %7.3f   %7.3f   %7.3f\n", Ra[i, k], Za[i, k], Ex[i, k])
            end
        end
    end

    #   Conditions to be satisfied by resulting equilibrium

    elxp = abs(Offset - Zxpsn) / Apl
    trixp = (Rmpl - Rxpsn) / Apl

    ityp = zeros(Int, 6)
    Rc = zeros(Nmax)
    Zc = zeros(Nmax)
    Rcc = zeros(16)
    Zcc = zeros(16)

    for j in 1:3
        ityp[j] = j
        Rc[j] = Rxpsn
        Zc[j] = Zxpsn
    end

    ityp[4] = 1
    Rc[4] = Rmpl - Apl
    Zc[4] = Offset
    ityp[5] = 1
    Rc[5] = Rmpl + Apl
    Zc[5] = Offset
    ityp[6] = 1
    Rc[6] = Rmpl - Apl * tri
    Zc[6] = Offset + Apl * El

    for j in 1:4
        ang = j * π / 10.
        Rcc[j] = Rmpl + Apl * cos(ang + tri * sin(ang))
        Zcc[j] = Offset + El * Apl * sin(ang)
        Rcc[j + 4] = Rmpl + Apl * cos(ang + π / 2. + tri * sin(ang + π / 2.))
        Zcc[j + 4] = Offset + El * Apl * sin(ang + π / 2.)
    end

    for j in 1:4
        al1 = Apl * (((1. + trixp) * (1. + trixp)) + elxp * elxp) / (2. * (1. + trixp))
        al2 = Apl * (((1. - trixp) * (1.0 - trixp)) + elxp * elxp) / (2. * (1.0 - trixp))
        anga = atan(2. * elxp * (1. + trixp) / (elxp * elxp - (1. + trixp) * (1. + trixp)))
        angb = atan(2. * elxp * (1. - trixp) / (elxp * elxp - (1. - trixp) * (1. - trixp)))
        rc1 = Rmpl + Apl - al1
        rc2 = Rmpl - Apl + al2
        ang1 = anga * j / 5.
        ang2 = angb * j / 5.
        Rcc[j + 8] = rc1 + al1 * cos(ang1)
        Zcc[j + 8] = Offset - al1 * sin(ang1)
        Rcc[j + 12] = rc2 - al2 * cos(ang2)
        Zcc[j + 12] = Offset - al2 * sin(ang2)
    end

    if (mprfg)
        println("  Single null case: boundary points ")
        for j in 1:16
            @printf(" %4i  %7.3f  %7.3f\n", j,Rcc[j],Zcc[j])
        end
    end

    expsi = zeros(Mr, Nz)
    fool = zeros(MN)
    psiext = zeros(MN, Mc)

    bb = zeros(Nmax, Mc)
    eb = zeros(llmax, Mc)
    cl = zeros(Mc + 1, Mc)
    fk = zeros(Mc + Nmax)

    for kk in 1:Mmax
        icl = ic[kk]
        for i in 1:Nz
            nof = (i - 1) * Mr
            for j in 1:Mr
                jn = nof + j
                expsi[jn] = 0.0
            end
        end

        for i in 1:icl
            if ( !(((Zmax - Za[i, kk]) * (Zmin - Za[i, kk]) <= 0.) && ((Rmax - Ra[i, kk]) * (Rmin - Ra[i, kk]) <= 0.)))
                for k = 1:Nm1:Nz
                    nof = (k - 1) * Mr
                    for j in 1:Mr
                        jn = nof + j
                        expsi[jn] = expsi[jn] + Ex[i, kk] * gfl(R[j], Ra[i, kk], Z[k] - Za[i, kk], 0.)
                    end
                end

                for k in 1:Nz
                    nof = (k - 1) * Mr
                    for j in 1:Mm1:Mr
                        jn = nof + j
                        expsi[jn] = expsi[jn] + Ex[i, kk] * gfl(R[j], Ra[i, kk], Z[k] - Za[i, kk], 0.)
                    end
                end
            end
        end

        eqsil!(expsi)

        for i in 1:icl
            if ( !(((Zmax - Za[i, kk]) * (Zmin - Za[i, kk]) > 0.) || ((Rmax - Ra[i, kk]) * (Rmin - Ra[i, kk]) > 0.)))
                for k in 1:Nz
                    nof = (k - 1) * Mr
                    for j in 1:Mr
                        jn = nof + j
                        expsi[jn] = expsi[jn] + Ex[i, kk] * gfl(R[j], Ra[i, kk], Z[k] - Za[i, kk])
                    end
                end
            end
        end

        for j in 1:MN           
            psiext[j, kk] = expsi[j]
        end

        # Computation of matrix elements for exact conditions

        splnco!(expsi)

        for j in 1:Nmax
            condit!(expsi, Rc[j], Zc[j], ityp[j], bb, j, kk)
        end
        for k in 1:llmax
            condit!(expsi, Rcc[k], Zcc[k], 1, eb, k, kk)
        end    
    
        # Computation of inductances

        cl[kk, kk] = 0.
        cl[Mmax + 1, kk] = 0.
        for i in 1:icl
            cl[Mmax + 1, kk] = cl[Mmax + 1, kk] + Ex[i, kk]
            cl[kk, kk] = cl[kk, kk] + Ex[i, kk]^2 * 1.0e6 * (0.58 + log(Ra[i, kk] / Rl[i, kk])) / (2. * π)
        end


        for i in 1:icl
            ii = i + 1
            if (ii <= ic[kk])
                for j in ii:icl
                    cl[kk, kk] = cl[kk, kk] + 2. * Ex[i, kk] * Ex[j, kk] * gfl(Ra[j, kk], Ra[i, kk], Za[j, kk] - Za[i, kk], 0.)
                end
            end
        end

        lp1 = kk + 1
        if (lp1 <= Mmax)
            for k in lp1:Mmax
                icm = ic[k]
                cl[kk, k] = 0.
                for i in 1:icl
                    for j in 1:icm
                        cl[kk, k] += Ex[i, kk] * Ex[j, k] * gfl(Ra[j, k], Ra[i, kk], Za[j, k] - Za[i, kk], 0.)
                    end
                end
            end
        end

    end
    
    # Computation of a new case

    icops, value = readdlm(IOBuffer(readline(dataFile)))
    icops = Int(icops)
    mpnmax = Mmax + Nmax + 1
    if (icops > 2)
        mpnmax += 1
    end
    totcurr, betapol, alfac = readdlm(IOBuffer(readline(dataFile)))
    raxis, zaxis, zdes, alp = readdlm(IOBuffer(readline(dataFile)))
    close(dataFile)

    jdes = floor(Int, 0.1 + (raxis - R[1]) / dr) + 1
    jaxis = jdes
    raxis = R[jaxis]
    naxis = floor(Int, 0.1 + (zaxis - Z[1]) / dz) + 1
    zaxis = Z[naxis]
    ndes = floor(Int, 0.1 + (abs(zdes) - Z[1]) / dz) + 1
    if (zdes > 0.)
        zdes = Z[ndes]
    end

    @printf("Magnetic Axis  r = %8.3f z = %8.3f\n", raxis,zaxis)
    @printf("Rail limiter   z = %8.3f alpfactor = %10.3e\n",zdes,alp)

    if (llmax > 0)
        alph = alp * 2. * π / (llmax * raxis)
    end

    close(dataFile)

    g = zeros(Mr, Nz)

    startt!()

    # Begin Iterations

    icycle = 1
    idecis = 0

    com = zeros(Mr)

    irsp = 0
    izsp = 0
    psicon = 0.
    fabs = 0.


    while (icycle <= 20)
        println(" ==== Cycle number ", icycle, " ====")

        eqsil!(g)

        compar!()

        for j in 1:MN
            fool[j] = 0.
            expsi[j] = g[j]
        end

        splnco!(expsi)

        for j in 1:Nmax
            condit!(expsi, Rc[j], Zc[j], ityp[j], bb, j, Mmax + 1)
        end
        for ll in 1:llmax
            condit!(expsi, Rcc[ll], Zcc[ll], 1, eb, ll, Mmax + 1)
        end

        xcur!()

        xt1 = 0.
        xt2 = 0.
        xt3 = 0.

        for i in 1:Mmax
            for j in 1:ic[i]
                curr = Ex[j,i] * fk[i]
                xt1 += curr^2
                xt2 += abs(curr)
                xt3 += abs(curr * Ra[j, i])
            end
        end

        if (mprfg)
            @printf(" SIG(I^2) = %12.4f SIG(ABS(I) = %12.4f  SIG(ABS(RI))) = %12.4f\n",xt1,xt2,xt3)
        end

        for k in 1:Mmax
            for j in 1:MN
                fool[j] += fk[k] * psiext[j,k]
                g[j] += fk[k] * psiext[j,k]
            end
        end

        saddle!()

        if (irsp > 2)
            if(mprfg)
                println(" Saddle point r = ", R[irsp], " z = ", Z[izsp])
            end
        end
        if (idecis > 0)
            break
        end

        curf!()

        global icycle += 1

    end

    plotit()

end    
