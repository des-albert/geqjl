function startt!()

    cp = 3. * totcurr / (20. * dr*dz) / abs(ndes - naxis)

    println(" Start current : ", totcurr)

    nmin = 2 * naxis - ndes
    njmax = ndes

    if(njmax < nmin)
        njmax = nmin
        nmin = ndes
    end
    for j in nmin:njmax
        for i in 1:5
            g[jaxis - 3 + i, j] = (cp * (1.0 - ((j - jaxis) / (ndes - naxis))^2)) * dz^2
        end
    end
    return nothing
end

