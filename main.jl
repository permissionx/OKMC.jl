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
            Record!(universe)
        end
        @do_every 1E7 quote
            print(universe)
            Dump(universe, dumpName)
        end
    end
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
const logName = "/mnt/c/Users/XUKE/Desktop/run.log"
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


