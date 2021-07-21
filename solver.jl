using Base: Float64
    using DelimitedFiles
    using SpecialFunctions

    include("bound_matrix.jl")
    include("eqsil.jl")
    include("flux.jl")

    Nexp = 6

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
        Zmax = Offset + 32 * (Apl*El) / 20
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

    for i in 1:Mr
        R[i] = Rmin + (i - 1) * dr
    end

    for j in 1:Nz
           Z[j] = Zmin + (j - 1) * dz
    end

    println("Rmin = ", Rmin, " Rmax = ", Rmax, " dr = ", dr)
    println("Zmin = ", Zmin, " Zmax = ", Zmax, " dz = ", dz)
    println("alpha = ", alpha, " sh = ", sh, " ss = ", ss)

    bound_matrix!()

    println("aux = ",aux[53, 12])


    # Read Poloidal Field Coil Data

    Mc = 15
    Ng = 1
    ic = zeros(Int, 15)
    Ra = zeros(Ng, Mc)
    Za = zeros(Ng, Mc)
    Ex = zeros(Ng, Mc)
    Rl = zeros(Ng, Mc)

    k = 0
    while (true)
        global k += 1
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

    close(dataFile)

    Mmax = k - 1

    if(mprfg)
        println("Conductor groups available for optimization")
        for k in 1:Mmax
            println("Group : ",k)
            for i in 1:ic[k]
                println("   ",Ra[i, k],"  ", Za[i, k],"   ",Ex[i, k])
            end
        end
    end

    #   Conditions to be satisfied by resulting equilibrium

    elxp = abs(Offset - Zxpsn)/Apl
    trixp = (Rmpl - Rxpsn)/Apl

    ityp = zeros(Int, 6)
    Rc = zeros(Mc)
    Zc = zeros(Mc)
    Rcc = zeros(16)
    Zcc = zeros(16)

    for j in 1:3
        ityp[j] = j
        Rc[j] = Rxpsn
        Zc[j] = Zxpsn
    end

    ityp[4] = 1
    Rc[4] = Rmpl - Apl
    Zc[4]= Offset
    ityp[5] = 1
    Rc[5] = Rmpl + Apl
    Zc[5] = Offset
    ityp[6] = 1
    Rc[6] = Rmpl - Apl*tri
    Zc[6]= Offset + Apl*El

    for j in 1:4
       ang = j*pi/10.
       Rcc[j] = Rmpl + Apl*cos(ang + tri*sin(ang))
       Zcc[j] = Offset + El*Apl*sin(ang)
       Rcc[j+4] = Rmpl + Apl*cos(ang + pi/2. + tri*sin(ang + pi/2.))
       Zcc[j+4] = Offset + El*Apl*sin(ang + pi/2.)
    end

    for j in 1:4
        al1 = Apl*(((1. + trixp)*(1. + trixp)) + elxp*elxp)/(2. * (1. + trixp))
        al2 = Apl*(((1. - trixp)*(1.0 - trixp)) + elxp*elxp)/(2. * (1.0 - trixp))
        anga = atan(2. * elxp*(1. + trixp)/(elxp*elxp - (1. + trixp)*(1. + trixp)))
        angb = atan(2. * elxp*(1. - trixp)/(elxp*elxp - (1. - trixp)*(1. - trixp)))
        rc1 = Rmpl + Apl - al1
        rc2 = Rmpl - Apl + al2
        ang1 = anga*j/5.
        ang2 = angb*j/5.
        Rcc[j+8] = rc1 + al1*cos(ang1)
        Zcc[j+8] = Offset - al1*sin(ang1)
        Rcc[j+12] = rc2 - al2*cos(ang2)
        Zcc[j+12] = Offset - al2*sin(ang2)
    end

    if (mprfg)
        println("  Single null case: boundary points ")
        for j in 1:16
            println(" ",j,"  ",Rcc[j],"  ",Zcc[j])
        end
    end

    expsi = zeros(MN)
    fool = zeros(MN)
    psiext = zeros(MN, Mc)

    for kk in 1:Mmax
        icl = ic[kk]
        for i in 1:Nz
           nof = (i - 1)*Mr
           for j in 1:Mr
              jn = nof + j
              expsi[jn] = 0.0
           end
        end

        for i in 1:icl
            if (! (((Zmax - Za[i,kk]))*(Zmin - Za[i,kk]) <= 0.) && ((Rmax - Ra[i,kk]))*(Rmin - Ra[i,kk]) <= 0.)
              for k = 1:Nm1:Nz
                 nof = (k - 1)*Mr
                 for j in 1:Mr
                    jn = nof + j
                    expsi[jn] = expsi[jn] + Ex[i,kk] * gfl(R[j], Ra[i,kk], Z[k] - Za[i,kk])
                 end
              end

              for k in 1:Nz
                 nof = (k - 1)*Mr
                 for j in 1:Mm1:Mr
                    jn = nof + j
                    expsi[jn] = expsi[jn] + Ex[i,kk] * gfl(R[j], Ra[i,kk], Z[k] - Za[i,kk])
                 end
              end
            end
        end

        eqsil!(expsi)
    end