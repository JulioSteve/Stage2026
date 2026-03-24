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

println("The total sum of probability = 1, we calculate the rest of the series after Ncut (from Ncut to infinity summing probabilities).\n")

for ncut in range(n0+1, n0+1+10)
    println("tested Ncut: ", ncut, " error/validity: ", round((P(n0, ncut))*100, digits=2)," / ",round((1.0-P(n0, ncut))*100, digits=2), "%")
end