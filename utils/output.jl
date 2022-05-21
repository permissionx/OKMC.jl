using UnicodePlots
#Dump uverserse in lammps dump format
function Dump(universe::Universe, fileName::String, mode::String="a")
    file = open(fileName, mode)
    write(file, "ITEM: TIMESTEP\n")
    write(file, "$(universe.nStep)\n")
    write(file, "ITEM: NUMBER OF ATOMS\n")
    write(file, "$(length(universe.defects))\n")
    write(file, "ITEM: BOX BOUNDS pp pp pp\n")
    write(file, "0 $(universe.mapSize[1])\n")
    write(file, "0 $(universe.mapSize[2])\n")
    write(file, "0 $(universe.mapSize[3])\n")
    write(file, "ITEM: ATOMS id type x y z size d1 d2 d3 radius\n")
    for defect in universe.defects
        if defect.type == 1
            direction = SIA_DIRECTIONS[defect.directionIndex]
        else
            direction = [0,0,0]
        end
        write(file, 
        "$(defect.index) $(defect.type) $(defect.coord[1]) $(defect.coord[2]) $(defect.coord[3]) \
        $(defect.size) $(direction[1]) $(direction[2]) $(direction[3]) $(defect.radius)\n")
    end
    close(file)
end

function RefreshFile!(fileName)
    file = open(fileName, "w")
    close(file)
end

macro print_distribution(p::Symbol)
    type = :(typeof(universe.defects[1].$p))
    name = String(p)
    return quote 
        ps = $(type)[]
        for defect in universe.defects
            push!(ps, defect.$p)
        end
        #println(ps)
        plot = histogram(ps, nbins=5)
        title!(plot, "Defect $($name)")
        xlabel!(plot, "")
        println(IOContext(stdout, :color=>true), plot)
    end
end



        
function RecordSV!(universe::Universe)
    nSia = 0
    nVac = 0
    for defect in universe.defects
        if defect.type == 1
            nSia += 1
        else
            nVac += 1
        end
    end
    push!(universe.history.nsSia, nSia)
    push!(universe.history.nsVac, nVac)
    push!(universe.history.steps, universe.nStep)
end


function PrintHistory(xs::Vector{Int64}, yses::Vector{Vector{T}}, names::Vector{String}; 
                      yscale=:identity) where T<:Real
    plot = lineplot(xs, yses[1], yscale=yscale, name = names[1], ylim = (0,maximum([maximum(yses[1]),maximum(yses[2])])+10))
    for i in 2:length(yses)
        lineplot!(plot, xs, yses[i], name=names[i])
    end
    title!(plot, "SIA/Vacancy number")
    xlabel!(plot, "Step")
    println(IOContext(stdout, :color=>true), plot)
end


function Base.print(symbol::Symbol, string::String)
    print(Crayon(foreground = symbol), string, Crayon(reset=true))
end


function Base.print(universe::Universe)
    #println("Step:", universe.nStep)
    run(`tput rc`)
    for _ in 1:OUTPUT_HEIGHTS
        println("                                                                        ") 
    end 
    run(`tput cuu $(OUTPUT_HEIGHTS)`)
    run(`tput sc`)
    print("🚀 ")
    print(:blue, "Step ")
    println(universe.nStep)
    print("👾 ")
    println(:red, "Defect number ")
    println("$(length(universe.defects)) including\
               \n $(universe.history.nsSia[end]) SIAs & $(universe.history.nsVac[end]) Vacancies)")
    nSia, nVac = SiaAndVacNumber(universe)
    print("$(nSia) single SIAs, $(nVac) single Vacancies")
    println("")
    print("📊 ")
    print(:green, "Current distributions\n")
    @print_distribution radius
    print("📅 ")
    print(:yellow, "Attributs\n")
    PrintHistory(universe.history.steps, [universe.history.nsVac, universe.history.nsSia], 
                ["Vacancy","SIA"])
    for _ in 1:2
        println()
    end
    flush(stdout)
end

function SiaAndVacNumber(universe::Universe)
    nSia = 0
    nVac = 0
    for defect in universe.defects
        if defect.type == 1
            nSia += defect.size
        else
            nVac += defect.size
        end
    end
    nSia, nVac
end
