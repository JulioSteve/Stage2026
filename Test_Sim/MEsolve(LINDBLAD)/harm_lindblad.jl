ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings, Printf, Base.Threads
theme(:dao)
palette = theme_palette(:dao)

const plot_lock = ReentrantLock()

function N_BE(x)
    if x == -1
        return 0  # χ→+∞ ≡ T→0 ⇒ N_th→0
    else
        return 1/(exp(x)-1)
    end
end

println("---------------------")
println("Threads disponibles : ", Threads.nthreads())
Prct = 0:1:100
function run_sim(Prct)
    n = length(Prct)
    Err = zeros(n)
    Time_iter = zeros(n)
        # Simulation Parameters
    # ===========================================
    α0 = 2.0 # Starting level for the initial coherent state
    χ = 1.0
    Nth = N_BE(χ) # mean thermal photon number in the cavity 
    tlist = range(0,10,1000)
    ω_k = 1.0 # unitless frequency system oscillation / coupling ratio (weak coupling regime if ω_k>>1)
    
    Threads.@threads for k in eachindex(Prct)
        Ncut = Nth + 4*√(Nth*(1+Nth)) + abs(α0)^2+4*abs(α0) # Dimension cutoff in the operators: μ+4σ of thermal state (steady state) + μ+4σ of initial coherent state

        prct = Prct[k]
        factor = 1.0+prct/100
        Ncut = Int(ceil(Ncut*factor)) # We take +prct% of the maximal population and round it to the upper value
        # println("⚠ Nth≈$(Int(ceil(Nth))) => Ncut=$Ncut")

        # Operators
        a = destroy(Ncut)

        c_decay = √(1+Nth)*a
        c_excite = √(Nth)*a'
        c_ops = [c_decay, c_excite]

        ψ0 = coherent(Ncut, α0)

        n_op = a'*a
        e_ops=[n_op]
        H = ω_k*a'*a 

        t_iter = @elapsed begin # calculate the time to solve
            sol = mesolve(H, ψ0, tlist, c_ops; e_ops=e_ops, progress_bar=Val(false)) # solver
        end
        n_anal(t) = α0^2*exp(-t)+Nth*(1-exp(-t)) # analytical solution
        N = n_anal.(tlist)
        Y = real.(sol.expect[1,:])
        Y2 = abs.(N.-Y)
        
        title = L"ħω=%$(χ)\times k_{B}T \qquad N_{cut}=%$(Ncut)"
        if  χ==-1
            title = L"T=0 \qquad N_{cut}=%$(Ncut)"
        end
        if prct in [0,25,50,75,100]
            lock(plot_lock) do
                p = plot(title=title*"(+$(prct)%)", legend=:topright)
                # Graph of the population
                plot!(p, tlist, Y,label="Simulation",lw=3, color=palette[1]);
                plot!(p, tlist, N, lw=2, label="Analytical model", color=palette[2], ls=:dash);
                xlabel!(p, L"$\tau=k t$ (unitless)")
                ylabel!(p, L"\langle \hat{n}(\tau)\rangle")

                # Error graph (inset inside the previous graph)
                yticks_val = range(minimum(Y2),maximum(Y2),5)
                yticks_label = [@sprintf("%.1e", v) for v in yticks_val]
                yticks = (yticks_val, yticks_label)

                plot!(p, tlist, Y2, color=palette[3], label="", yticks=yticks, title="Absolute error", lw=2 ; inset=bbox(0.6, 0.4, 0.4, 0.4), subplot=2)
                savefig(p, "$(prct)%.svg")
            end
        end
        Err[k] = maximum(Y2)
        Time_iter[k] = t_iter
    end
    println()
    return Err, Time_iter
end

function lastgraph(Prct, Err, T)
    p = plot(title="Maximal error VS cutoff offset(%)", right_margin=-50Plots.mm)
    
    yticks = 10.0 .^ (floor(Int, log10(minimum(Err))):ceil(Int, log10(maximum(Err))))
    plot!(p, Prct, Err, label="", yticks=yticks, yscale=:log10, color=palette[1], lw=3, xlabel="Truncation offset (+%)", ylabel="Max absolute error") 

    T .*= 1000 # conversion en ms
    T[1] = minimum(T)
    yticks = round.(range(minimum(T), maximum(T), 5), digits=0)
    
    plot!(p, Prct, T, color=palette[2], label="", yticks=Int.(yticks), title="Computation time", lw=2, ylabel="Time (ms)"; inset=bbox(0.5, 0.15, 0.4, 0.4), subplot=2)
    
    savefig(p, "stability.svg")
end

Errors, Times = run_sim(Prct)
println("---------------------")
lastgraph(Prct, Errors, Times)