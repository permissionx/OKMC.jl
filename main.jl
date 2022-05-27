include("utils/kmc.jl")

function Run!(universe::Universe)
    Init!(universe)
    while universe.nStep <= 1E11
        if universe.nStep % 1E4 == 0
            defect = Defect(rand(0.:199.,3), rand(1:2), rand(1:4), rand(1:1))
            push!(universe, defect)
        end
        #exit()
        IterStep!(universe)
        if universe.nStep % 1E5 == 0
            Record!(universe)
        end
        if universe.nStep % 1E7 == 0
            print(universe)
            Dump(universe, dumpName)
        end
    end
end

# -----------------------INPUT---------------------------
Random.seed!(31415926)
const SIA_DISAPPEAR_RATE = 1.4E-7
const VAC_DISAPPEAR_RATE = 4.5E-8
const MAX_DEFECT_SIZE = 30000
const OUTPUT_HEIGHTS = 40
mapSize = [200.,200.,200.]
cellLength = 20
const k_B = 8.617E-5 # eV/K
const VAC_MIGRATE_BARRIERS = [1.66, 1.66, 0.89, 1.56]  # eV
const VAC_K0 = 6E12 # s
const VAC_BINDING_ENERGYS = [0.0, -0.12, 0.01, 0.64, 0.52, 1.21, 0.77, 0.96, 0.85, 1.56, 0.96, 1.56, 1.21, 2.49, 2.62,
                             0.96, 1.00, 0.85, 0.88, 1.59, 1.76, 1.79, 2.14, 2.56]
const SIA_K0 = 1.43E12
const SIA_ALPHA = -0.365
const SIA_MIGRATE_BARRIER = 0.04
const SIA_STEER_BARRIER = 0.38
const TEMPERATURE = 1000.0 # K
const dumpName = "./run.dump"
const logName = "./run.log"
# -----------------------END INPUT-------------------------

universe = Universe(mapSize, cellLength)
#Run!(universe)
#Run!(universe)
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
# realistic frequency ✔️  check carefully ❓
# beatifify screen output ✔️
# outpot dataframe for python plot ✔️

