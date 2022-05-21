include("utils/kmc.jl")

function Run!(universe::Universe)
    Init!(universe)
    while universe.nStep <= 1E11
        @do_every 1E4 quote
            defect = Defect(rand(0.:199.,3), 2, rand(1:4), rand(1:1))
            push!(universe, defect)
        end
        #exit()
        IterStep!(universe)
        @do_every 1E6 RecordSV!(universe)
        @do_every 1E7 quote
            print(universe)
            Dump(universe, dumpName)
        end
    end
end



Random.seed!(31415926)
const SIA_DISAPPEAR_RATE = 0
const MAX_DEFECT_SIZE = 1000
const OUTPUT_HEIGHTS = 40
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


