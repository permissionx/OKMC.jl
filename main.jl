include("utils/kmc.jl")

function Run!(universe::Universe)
    Init!(universe)
    while universe.nStep <= 1E10
        @do_every 1E4 quote
            defect = Defect(rand(0.:199.,3), rand(1:2), rand(1:4), rand(1:1))
            push!(universe, defect)
        end
        #exit()
        IterStep!(universe)
        @do_every 1E5 quote 
            RecordSV!(universe)
            Record!(universe)
        end
        @do_every 1E7 quote
            print(universe)
            Dump(universe, dumpName)
        end
    end
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
    print(:red, "Defect number $(length(universe.defects))")
    print(" including\n")
    println(" $(universe.history.nsSia[end]) SIAs & $(universe.history.nsVac[end]) Vacancies")
    nSia, nVac = SiaAndVacNumber(universe)
    print(" $(nSia) single SIAs & $(nVac) single Vacancies")
    println("")
    print("📊 ")
    print(:green, "Current distributions\n")
    @print_distribution radius
    print("📅 ")
    print(:yellow, "Attributs\n")
    PrintHistory(universe.history.steps, [universe.history.nsVac, universe.history.nsSia], 
                ["Vacancy","SIA"], "SIA/vacancy number")
    #PrintHistory(universe.history.steps, [universe.history.sinkedSiasVector, universe.history.annihilatedSiasVector], 
    #            ["sinked","annihilated"], "Sinked/annihilated SIA")
    for _ in 1:2
        println()
    end
    flush(stdout)
end

Random.seed!(31415926)
const SIA_DISAPPEAR_RATE = 0.0001
const MAX_DEFECT_SIZE = 10000
const OUTPUT_HEIGHTS = 50
mapSize = [200.,200.,200.]
cellLength = 20
universe = Universe(mapSize, cellLength)
#const dumpName = "./run/run.dump"
const dumpName = "/mnt/c/Users/XUKE/Desktop/run.dump"
RefreshFile!(dumpName)
Run!(universe)
#InitRadius!(universe)
#Run!(universe)
#Init!(universe)
#Init!(universe)
#Run!(universe)
#empty
#using Profile, PProf
#@profile Run!(universe)
#pprof(;webport=58599)


# todo: 
# fix boundary cells ✔️
# realistic probability ❓
# beatifify screen output ✔️
# outpot dataframe for python plot ✔️


