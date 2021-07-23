function saddle!()

    f1(a,b) = (a - hk)*(b - hk)

    irsp = 0
    for i in 2:Nm1
        for k in 2::Nm1
            hk = g[k,i]
            if ( !(irsp != 0 && psicon >= hk))
                zfrr = g[k + 1, i] - 2. * hk + g[k - 1,i]
                zfzz = g[k, i+1] - 2. * hk + g[k, i - 1]
                zfrz = g[k + 1, i + 1] + g[k - 1, i - 1] - g[k + 1, i - 1] - g[k - 1, i + 1]
                if (  (16.*zfrr * zfrr * zfzz - zfrz^2) < 0. )
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
                        if ( f1(g[k + ko, i + io], g[k - ko, i - io]) > 0.)
                            if ( g[k + ko, i + io] != hk)
                                sg = true
                            else
                                gg = true
                            end
                            if (sg && gg)
                                @goto finish
                            end
                        end
                    end
                    @goto last
@label finish
                    irsp = k
                    izsp = i
                    psicon = hk
                end
@label last
            end
        end
    end

    return nothing

end