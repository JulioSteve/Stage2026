n0 = 5.0

function P(nmoy, Ncut)
    eps = 0.0
    pn = exp(-nmoy) * nmoy^Ncut / factorial(big(Int(Ncut)))
    for n in Ncut:Int(1e3)
        eps += pn
        pn *= nmoy / (n + 1)
    end
    
    return round(Float64(eps), digits=5)
end

P(n0, 7)