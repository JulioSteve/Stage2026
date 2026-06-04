ENV["GKSwstype"] = "100"
using QuantumToolbox, Plots, LaTeXStrings, Printf, Base.Threads
theme(:dao)
palette = theme_palette(:dao)

#######################

χ=Inf

#######################

function BoseEinstein(χ)
    """
    Distribution de Bose-Einstein avec la notation χ=ħω/(k_B T).
    Le cas T->0 est équivalent à χ=+∞ et donc il faut évaluer la fonction en écrivant "Inf" (géré par Julia proprement).
    """
    return 1.0/(exp(χ)-1.0)
end

function which_system(name::String)
    Nth = BoseEinstein(χ)

    if name=="Qubit"
        
    elseif name=="Harmonic"
        println("Harmonic WIP")
    else
        println("Write either 'Qubit' or 'Harmonic'.")
        return Nothing
    end
end


####################################



####################################

