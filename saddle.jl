function saddle!()

    global psicon, irsp, izsp

    f1(a,b,c) = (a - c)*(b - c)
   
    irsp = 0
    for i in 2:Nm1
        for k in 2:Mm1
            hk = g[k,i]
            if ( !(irsp != 0 && psicon >= hk))
                zfrr = g[k + 1, i] - 2. * hk + g[k - 1,i]
                zfzz = g[k, i+1] - 2. * hk + g[k, i - 1]
                zfrz = g[k + 1, i + 1] + g[k - 1, i - 1] - g[k + 1, i - 1] - g[k - 1, i + 1]
                if (  (16. * zfrr * zfzz - zfrz^2) < 0. )
                    sg = false
                    gg = false
                    ko = 1
                    io = -2
                    for l in 1:4
                        if (l < 4)
                            io += 1
                        end
                        if (l == 4)
                            ko = 0
                        end
                        if ( f1(g[k + ko, i + io], g[k - ko, i - io], hk) > 0.)
                            if ( g[k + ko, i + io] > hk)
                                sg = true
                            else
                                gg = true
                            end
                            if (sg && gg)
                                @goto L10
                            end
                        end
                    end
                    @goto L20
@label L10
                    irsp = k
                    izsp = i
                    psicon = hk
                end
            end
@label L20
        end
    end

    return nothing

end