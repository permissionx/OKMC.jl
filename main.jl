include("head.jl")
include("output.jl")

function RefreshFile!(fileName)
    file = open(fileName, "w")
    close(file)
end

CellCoord(universe::Universe, coord::Vector{Int64}) = floor.(Int64, coord/universe.cellLength ) .+ 1

function GetCell(universe::Universe, cellCoord::Vector{Int64})
    for i in 1:3
        if cellCoord[i] > universe.nsCells[i]
            cellCoord[i] -= 1
        end
    end
    return universe.cells[cellCoord[1], cellCoord[2], cellCoord[3]]
end

function PBCCellCoord!(universe::Universe, cellCoord::Vector{Int64})
    crossSign = Int64[]
    for i in 1:3
        if cellCoord[i] < 1
            cellCoord[i] = universe.nsCells[i] + cellCoord[i]
            push!(crossSign, -1)
        elseif cellCoord[i] > universe.nsCells[i]
            cellCoord[i] = cellCoord[i] - universe.nsCells[i]
            push!(crossSign, 1)
        else
            push!(crossSign, 0)
        end
    end
    return crossSign
end

function InitCells!(universe::Universe)
    for i in 1:universe.nsCells[1]
        for j in 1:universe.nsCells[2]
            for k in 1:universe.nsCells[3]
                cell = universe.cells[i, j, k]
                for di in -1:1
                    for dj in -1:1
                        for dk in -1:1
                            neighborCellCoord = [i, j, k] + [di, dj, dk]
                            crossSign = PBCCellCoord!(universe, neighborCellCoord)
                            neighborCell = GetCell(universe, neighborCellCoord)
                            if crossSign == [0,0,0]
                                push!(cell.normalNeighbors, neighborCell)
                            else
                                push!(cell.crossNeighbors, neighborCell)
                                push!(cell.crossSigns, crossSign)
                            end
                        end
                    end
                end
            end
        end
    end
end

function Delta(universe::Universe, coord1::Vector{Int64}, coord2::Vector{Int64}) 
    # vactor points from 1 to 2 
    (coord2 - coord1)
end

function Delta(universe::Universe, coord1::Vector{Int64}, coord2::Vector{Int64}, crossSign::Vector{Int64})
    delta = coord2 - coord1
    for i in 1:3
        delta[i] += crossSign[i] * universe.mapSize[i]
    end
    delta
end

function Distance(universe::Universe, coord1::Vector{Int64}, coord2::Vector{Int64})
    delta = Delta(universe, coord1, coord2)
    sqrt(sum(delta .* delta))
end

function Distance(universe::Universe, coord1::Vector{Int64}, coord2::Vector{Int64}, crossSign::Vector{Int64})
    delta = Delta(universe, coord1, coord2, crossSign)
    sqrt(sum(delta .* delta))
end

function FindNeighbors(universe::Universe, defect::Defect)
    cell = universe.cells[defect.cellIndex]
    neighbors = Defect[]
    for cell in cell.normalNeighbors
        for neighbor in cell.defects
            if neighbor.index != defect.index
                if Distance(universe, defect.coord, neighbor.coord) <= neighbor.radius+defect.radius
                    push!(neighbors, neighbor)
                end
            end
        end
    end
    for (cell, crossSign) in zip(cell.crossNeighbors, cell.crossSigns)
        for neighbor in cell.defects
            if neighbor.index != defect.index
                if Distance(universe, defect.coord, neighbor.coord, crossSign) <= neighbor.radius+defect.radius
                    push!(neighbors, neighbor)
                end
            end
        end
    end
end

function PBCCoord!(universe::Universe, coord::Vector{Int64})
    for i in 1:3
        if coord[i] < 0
            coord[i] = coord[i] + universe.mapSize[i]
        elseif coord[i] >= universe.mapSize[i]
            coord[i] = coord[i] - universe.mapSize[i]
        end
    end
end

function Reaction!(universe::Universe, defect1::Defect, defect2::Defect)
    largeDefect = (defect1.size > defect2.size) ? defect1 : defect2
    smallDefect = (defect1.size <= defect2.size) ? defect1 : defect2
    if largeDefect.type == smallDefect.type
        #combine
        center = (largeDefect.coord*largeDefect.size + smallDefect.coord*smallDefect.size) / (largeDefect.size + smallDefect.size)
        newCoord = round.(Int64, center)
        ChangeSize!(universe, largeDefect, largeDefect.size + smallDefect.size)
        delete!(universe, smallDefect)
        Move!(universe, largeDefect, newCoord)
    else
        if largeDefect.size == smallDefect.size
            delete!(universe, smallDefect)
            delete!(universe, largeDefect)
        else
        #swallow
            radius = largeDefect.radius 
            ChangeSize!(universe, largeDefect, largeDefect.size - smallDefect.size)
            moveLength = radius - largeDefect.radius
            delta = Delta(universe, largeDefect.coord, smallDefect.coord)
            distance = sqrt(sum(delta .* delta))
            delete!(universe, smallDefect)
            if distance != 0
                newCoordFloat = Vector{Float64}(undef,3)
                for i in 1:3
                    newCoordFloat[i] = largeDefect.coord[i] - moveLength * delta[i] / distance
                end
                newCoord = round.(Int64, newCoordFloat)
                Move!(universe, largeDefect, newCoord)
            end
        end
    end
end


function Changed!(universe::Universe, defect::Defect)
    neighbors = FindNeighbors(universe, defect)
    if length(neighbors) > 0 
        neighbor = sample(neighbors)
        Reaction!(universe, defect, neighbor)
    end
end


function SiaProbabilities(size::Int64)
    # 1 for migration, 2 for sterring
    return [1.,0.0]
end


function VacProbabilities(size::Int64)
    # 1 for migration, 2 for emittion
    if size == 1
        return Float64[0.5,0.0]
    else
        return Float64[0,0.0000]
    end
end


function ChangeSize!(universe::Universe, defect::Defect, size::Int64)
    defect.size = size
    Radius!(universe, defect)
    ChangeBehaviors!(universe, defect)
end

function Radius!(universe::Universe, defect::Defect)
    if defect.type == 1
        defect.radius = universe.constants.siaRadii[defect.size]
    else
        defect.radius = universe.constants.vacRadii[defect.size]
    end
end

function InitBehaviors!(universe::Universe, defect::Defect)
    if defect.type == 1
        types = [1, 3]
        probabilities = SiaProbabilities(defect.size)
    else
        types = [1, 2]
        probabilities = VacProbabilities(defect.size)
    end
    behaviors = Behaviors()
    behaviors.types = types
    behaviors.probabilities = probabilities
    probability = sum(probabilities)
    behaviors.probability = probability
    defect.behaviors = behaviors
    push!(universe.defectProbabilities, probability)
    universe.totalProbability += probability
end

function ChangeBehaviors!(universe::Universe, defect::Defect)
    oldProbability = defect.behaviors.probability
    if defect.type == 1
        probabilities = SiaProbabilities(defect.size)
    else
        probabilities = VacProbabilities(defect.size)
    end
    defect.behaviors.probabilities = probabilities
    probability = sum(probabilities)
    i = findfirst(x->x.index===defect.index, universe.defects)
    universe.defectProbabilities[i] = probability
    universe.totalProbability += probability - oldProbability
end



function Base.push!(universe::Universe, defect::Defect)
    cellCoord = CellCoord(universe, defect.coord)
    cell = GetCell(universe, cellCoord)
    push!(cell, defect)
    universe.maxIndex += 1
    defect.index = universe.maxIndex
    Radius!(universe, defect)
    InitBehaviors!(universe, defect)
    push!(universe.defects, defect)
    Changed!(universe, defect)
end


function Base.push!(cell::Cell, defect::Defect)
    push!(cell.defects, defect)
    defect.cellIndex = cell.index
end


function Base.delete!(universe::Universe, defect::Defect)
    cell = universe.cells[defect.cellIndex]
    delete!(cell, defect)
    i = findfirst(x->x.index==defect.index, universe.defects)
    deleteat!(universe.defects, i)
    deleteat!(universe.defectProbabilities, i)
    universe.totalProbability -= defect.behaviors.probability
end


function Base.delete!(cell::Cell, defect::Defect)
    deleteat!(cell.defects, findfirst(x->x.index==defect.index, cell.defects))
end

function Move!(universe::Universe, defect::Defect, coord::Vector{Int64})
    PBCCoord!(universe, coord)
    defect.coord = coord
    CrossCells!(universe, defect)
    Changed!(universe, defect)
end

function CrossCells!(universe::Universe, defect::Defect)
    cellCoord = CellCoord(universe, defect.coord)
    cell = GetCell(universe, cellCoord)
    if cell.index != defect.cellIndex
        delete!(universe.cells[defect.cellIndex], defect)
        push!(cell, defect)
    end
end


#KMC
function RandomADefect(universe::Universe)
    weights = Weights(universe.defectProbabilities)
    sample(universe.defects, weights)
end

function RandomABehavior(defect::Defect)
    sample(defect.behaviors.types, Weights(defect.behaviors.probabilities))
end

function Behave!(universe::Universe, defect::Defect, type::Int64)
    if type == 1
        Migrate!(universe, defect)
    elseif type == 2
        Emit!(universe, defect)
    else
        Steer!(defect)
    end
end

function Migrate!(universe::Universe, defect::Defect)
    if defect.type == 1
        r = rand()
        if r <= SIA_DISAPPEAR_RATE
            delete!(universe, defect)
            return
        end
        sign = sample([1,-1])
        displace = SIA_DIRECTIONS[defect.directionIndex] * sign
        newCoord = defect.coord + displace
        Move!(universe, defect, newCoord)
    else
        displace = sample(FIRST_NEIGHBORS)
        newCoord = defect.coord + displace
        Move!(universe, defect, newCoord)
    end
end

function MigrateSecondNeighbor!(universe::Universe, defect::Defect)
    displace = sample(SECOND_NEIGHBORS)
    newCoord = defect.coord + displace
    Move!(universe, defect, newCoord)
end

function Emit!(universe::Universe, defect::Defect)
    r = defect.radius+2
    coord = RandomAPointOnShpereSurface(r) + defect.coord
    PBCCoord!(universe, coord)
    vac = Defect(coord, 2, 1, 1)
    ChangeSize!(universe, defect, defect.size-1)
    push!(universe, vac)    
end

function RandomAPointOnShpereSurface(r::Float64)
    # pick a random point on a sphere surface uniformly
    # http://mathworld.wolfram.com/SpherePointPicking.html
    u = rand()*2-1
    theta = 2pi*rand()
    r2d = sqrt(1-u*u)
    x = r2d * cos(theta) * r
    y = r2d * sin(theta) * r
    z = u*r
    return round.(Int64,[x,y,z])
end


function Steer!(defect::Defect)
    @assert(defect.type==1, "steering a vacancy!")
    newDirections = [a for a in Int64[1,2,3,4] if a != defect.directionIndex]
    newDirection = sample(newDirections)
    defect.directionIndex = newDirection
end


function IterStep!(universe::Universe)
    universe.nStep += 1
    defect = RandomADefect(universe)
    type = RandomABehavior(defect)
    Behave!(universe, defect, type)
end


function InitRadius!(universe::Universe)
    vacRadii = Float64[]
    for n in 1:MAX_DEFECT_SIZE
        push!(vacRadii, (3/4/pi*2^3/2*n)^(1/3) - (3/4/pi*2^3/2)^(1/3) + 3^(1/3))
    end
    siaRadii = Float64[]
    for n in 1:MAX_DEFECT_SIZE
        push!(siaRadii, (2^2/3^(1/2)/pi*n)^(1/2) - (2^2/3^(1/2)/pi)^(1/2) + 3^(1/3))
    end
    universe.constants.vacRadii = vacRadii
    universe.constants.siaRadii = siaRadii
end

function Init!(universe::Universe)
    InitCells!(universe)
    InitRadius!(universe)
    run(`tput sc`)
end

function End(universe::Universe)
    run(`tput sc`)
end

function Test!(universe::Universe)
    Init!(universe)
    defect = Defect(rand(0:99,3), 2, rand(1:4), rand(10:20))
    push!(universe, defect)
    Dump(universe, dumpName)
end

macro do_every(n::Int64, f::Expr)
    return quote 
        if universe.nStep%$n == 0
            eval($f)
        end
    end
end


function Run_small!(universe::Universe)
    Init!(universe)
    while universe.nStep <= 1_000
        @do_every 1 quote
            defect = Defect(rand(0:199,3), rand(1:2), rand(1:4), rand(1:1))
            push!(universe, defect)
        end
        #exit()
        IterStep!(universe)
        #@do_every 1 RecordSV!(universe)
        #@do_every 1 quote
        #    print(universe)
        #    Dump(universe, dumpName)
        #end
    end
end

function Run!(universe::Universe)
    Init!(universe)
    while universe.nStep <= 1_000_000_000_000
        @do_every 100_000 quote
            defect = Defect(rand(0:199,3), rand(1:2), rand(1:4), rand(1:1))
            push!(universe, defect)
        end
        #exit()
        IterStep!(universe)
        @do_every 10_000_000 RecordSV!(universe)
        @do_every 100_000_000 quote
            print(universe)
            Dump(universe, dumpName)
        end
    end
end



Random.seed!(31415926)
const SIA_DISAPPEAR_RATE = 1E-7
const MAX_DEFECT_SIZE = 10000
const OUTPUT_HEIGHTS = 40
mapSize = [200,200,200]
cellLength = 20
universe = Universe(mapSize, cellLength)
const dumpName = "./run/run.dump"
RefreshFile!(dumpName)
Run_small!(universe)

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
# fix boundary cells ❓
# realistic probability ❓
# beatifify screen output ✔️
# outpot dataframe for python plot ✔️


