function axslen(i1::Int64, i2::Int64)
ccall((:axslen, "libdislin_d.so"), Nothing, (Int32,Int32), i1, i2)
end
function color(s::String)
ccall((:color, "libdislin_d.so"), Nothing, (Ptr{UInt8},), s)
end
function complx()
ccall((:complx, "libdislin_d.so"), Nothing, ())
end
function contur(xray::Vector{Float64}, n::Int64, yray::Vector{Float64},
    m::Int64, zmat::Array{Float64, 2}, zlev::Float64)
    zt = transpose(zmat)
ccall((:contur, "libdislin_d.so"), Nothing, (Ptr{Float64}, Int32, Ptr{Float64}, 
    Int32, Ptr{Float64}, Float64), xray, n, yray, m, zt, zlev)
end
function disini()
ccall((:disini, "libdislin_d.so"), Nothing, () )
end
function disfin()
ccall((:disfin, "libdislin_d.so"), Nothing, ())
end
function endgrf()
ccall((:endgrf, "libdislin_d.so"), Nothing, ())
end
function graf(x1::Float64, x2::Float64, x3::Float64, x4::Float64,
    x5::Float64, x6::Float64, x7::Float64, x8::Float64)
ccall((:graf, "libdislin_d.so"), Nothing, (Float64, Float64, Float64,
    Float64, Float64, Float64, Float64, Float64), 
    x1, x2, x3, x4, x5, x6, x7, x8)
end
function pagera()
ccall((:pagera, "libdislin_d.so"), Nothing, ())
end
function pagfll(i::Int64)
ccall((:pagfll, "libdislin_d.so"), Nothing, (Int32,), i)
end
function rlrec(x1::Float64, x2::Float64, x3::Float64, x4::Float64)
ccall((:rlrec, "libdislin_d.so"), Nothing, (Float64, Float64, Float64, Float64), 
            x1, x2, x3, x4)
end
function setpag(s::String)
ccall((:setpag, "libdislin_d.so"), Nothing, (Ptr{UInt8},), s)
end
function xaxgit()
ccall((:xaxgit, "libdislin_d.so"), Nothing, ())
end

function plotit()

global psicon
global fabs


    setpag("ps4l")
    disini()
    complx()
    pagera()
    pagfll(255)
    color("blue")
    axslen(1120, 1600)
    graf(0., 14., 0., 2., -10., 10., -10., 2.)
    xaxgit()
    color("green")
    rlrec(1.3135,6.0605,.749,1.979)
    rlrec(1.3135,4.0325,.749,1.979)
    rlrec(1.3135,2.0035,.749,1.979)
    rlrec(1.3135,-0.0245,.749,1.979)
    rlrec(1.3135,-2.0535,.749,1.979)
    rlrec(1.3135,-4.0815,.749,1.979)
    rlrec(3.47,8.045,.968,.976)
    rlrec(7.9845,6.8275,.649,.595)
    rlrec(11.581,3.8275,.708,1.125)
    rlrec(11.5805,-1.6805,.649,1.125)
    rlrec(7.985,-6.2575,.82,.945)
    rlrec(3.4705,-7.069,1.633,.976)

    color("blue")

    for i in 2:16
        plev = (psicon + (i - 1)*(fabs - psicon)/15.)
        contur(R, Mr, Z, Nz, g, plev)
    end 
    color("red")
    contur(R, Mr, Z, Nz, g, psicon)
    endgrf()
    disfin()

end