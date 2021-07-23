function compar!()

    if (icycle >= 1)
        rel = 0.
        for j in 1:Mr
            tot = abs(g[j, naxis]) / 2. + abs(com[j]) / 2.
            dev = abs(g[j, naxis] - com[j])
            ren = dev / tot
            if (ren > rel)
                rel = ren
            end
        end
        println(" Relative Error = ", rel)
        if(rel < error)
            idecis = 1
            return nothing
        end
    end
    for j in 1:Mr
        com[j] = g[j, naxis]
    end
    idecis = 0
    return nothing
end