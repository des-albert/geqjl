function eqsil!(q)

    qq = copy(q)

    npn = Nm1*Mr
    for i in 1:Mr
        qq[i] = 0.
        qq[i + npn] = 0.
    end

    for k in 1:Mr:MN
        qq[k] = 0.
        qq[k + Mm1] = 0.
    end

    flux!(qq)

    for i in 1:llp
         sum = 0.
         for l = 1:llp
            sum += qq[ip[l]] * aux[l, i]
         end
         q[jp[i]] = q[jp[i]] + sum
    end

    flux!(q)

    return nothing

end
