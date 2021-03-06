include("head.jl")
include("output.jl")

CellCoord(universe::Universe, coord::Vector{Float64}) = floor.(Int64, coord/universe.cellLength ) .+ 1

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

function Delta(universe::Universe, coord1::Vector{Float64}, coord2::Vector{Float64}) 
    # vactor points from 1 to 2 
    (coord2 - coord1)
end

function Delta(universe::Universe, coord1::Vector{Float64}, coord2::Vector{Float64}, crossSign::Vector{Int64})
    delta = coord2 - coord1
    for i in 1:3
        delta[i] += crossSign[i] * universe.mapSize[i]
    end
    delta
end

function Distance(universe::Universe, coord1::Vector{Float64}, coord2::Vector{Float64})
    delta = Delta(universe, coord1, coord2)
    sqrt(sum(delta .* delta))
end

function Distance(universe::Universe, coord1::Vector{Float64}, coord2::Vector{Float64}, crossSign::Vector{Int64})
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
    crossSigns = [[0,0,0] for _ in 1:length(neighbors)]
    for (cell, crossSign) in zip(cell.crossNeighbors, cell.crossSigns)
        for neighbor in cell.defects
            if neighbor.index != defect.index
                if Distance(universe, defect.coord, neighbor.coord, crossSign) <= neighbor.radius+defect.radius
                    push!(neighbors, neighbor)
                    push!(crossSigns, crossSign)
                end
            end
        end
    end
    neighbors, crossSigns
end

function PBCCoord!(universe::Universe, coord::Vector{Float64})
    for i in 1:3
        if coord[i] < 0
            coord[i] = coord[i] + universe.mapSize[i]
        elseif coord[i] >= universe.mapSize[i]
            coord[i] = coord[i] - universe.mapSize[i]
        end
    end
end

function Reaction!(universe::Universe, defect1::Defect, defect2::Defect, crossSign::Vector{Int64})
    #if defect1.type == defect2.type == 1   # for no reaction for sia and sia
    #    return
    #end
    if defect1.size > defect2.size
        largeDefect = defect1
        largeDefectCoord = defect1.coord
        smallDefect = defect2
        smallDefectCoord = defect2.coord + crossSign.*universe.mapSize
    else
        largeDefect = defect2
        largeDefectCoord = defect2.coord+ crossSign.*universe.mapSize
        smallDefect = defect1
        smallDefectCoord = defect1.coord 
    end
    if largeDefect.type == smallDefect.type
        #combine
        newCoord = (largeDefectCoord*largeDefect.size + smallDefectCoord*smallDefect.size) / (largeDefect.size + smallDefect.size)
        ChangeSize!(universe, largeDefect, largeDefect.size + smallDefect.size)
        delete!(universe, smallDefect)
        Move!(universe, largeDefect, newCoord)
    else
        if largeDefect.size == smallDefect.size
            delete!(universe, smallDefect)
            delete!(universe, largeDefect)
        else
        #swallow: need to test on the boundary
            radius = largeDefect.radius 
            ChangeSize!(universe, largeDefect, largeDefect.size - smallDefect.size)
            moveLength = radius - largeDefect.radius
            delta = Delta(universe, largeDefectCoord, smallDefectCoord)
            distance = sqrt(sum(delta .* delta))
            delete!(universe, smallDefect)
            if distance != 0
                newCoord = Vector{Float64}(undef,3)
                for i in 1:3
                    newCoord[i] = largeDefect.coord[i] - moveLength * delta[i] / distance
                end
                Move!(universe, largeDefect, newCoord)
            end
        end
        universe.record.annihilatedSiaNum += smallDefect.size
    end
end


function Changed!(universe::Universe, defect::Defect)
    neighbors, crossSigns = FindNeighbors(universe, defect)
    if length(neighbors) > 0 
        i = rand(1:length(neighbors))
        neighbor = neighbors[i]
        crossSign = crossSigns[i]
        Reaction!(universe, defect, neighbor, crossSign)
    end
end


function SiaFrequencies(universe::Universe, size::Int64)
    # 1 for migration, 2 for sterring
    if size == 1
        return [universe.constants.siaMigrateFrequencies[size], universe.constants.siaSteerFrequency]
    else
        return [universe.constants.siaMigrateFrequencies[size], 0.0]
    end
end

function VacFrequencies(universe::Universe, size::Int64)
    # 1 for migration, 2 for emittion
    if size <= length(VAC_MIGRATE_BARRIERS)
        return [universe.constants.vacMigrateFrequencies[size], universe.constants.vacEmitFrequencies[size]]
    else
        return [0.0, universe.constants.vacEmitFrequencies[size]]
    end
end

function InitFrequencies!(universe::Universe)
    InitVacFrequencies!(universe::Universe)
    InitSiaFrequencies!(universe::Universe)
end

function InitSiaFrequencies!(universe::Universe)
    # units : s, eV
    universe.constants.siaMigrateFrequencies = Float64[]
    for i in 1:MAX_DEFECT_SIZE
        k0 = SIA_K0*i^SIA_ALPHA
        frequency = k0*exp(-SIA_MIGRATE_BARRIER/universe.temperature/k_B)
        push!(universe.constants.siaMigrateFrequencies, frequency)
    end
    # todo: check the k0 for the steering ???
    universe.constants.siaSteerFrequency = SIA_K0*exp(-SIA_MIGRATE_BARRIER/universe.temperature/k_B)
end



function InitVacFrequencies!(universe::Universe)
    # units : s, eV
    universe.constants.vacMigrateFrequencies = Float64[]
    universe.constants.vacEmitFrequencies = Float64[]
    k0 = VAC_K0
    for i in 1:MAX_DEFECT_SIZE
        if i <= length(VAC_MIGRATE_BARRIERS)
            barrier = VAC_MIGRATE_BARRIERS[i]
            frequency = k0*exp(-barrier/universe.temperature/k_B)
            push!(universe.constants.vacMigrateFrequencies, frequency)
        else
            push!(universe.constants.vacMigrateFrequencies, 0.0)
        end
    end

    # todo: check the k0 for the emit ???
    for i in 1:MAX_DEFECT_SIZE
        if i == 1
            push!(universe.constants.vacEmitFrequencies, 0.0)
        else
            if i<= length(VAC_BINDING_ENERGYS)
                barrier = VAC_BINDING_ENERGYS[i] + VAC_MIGRATE_BARRIERS[1]
            else
                barrier = VAC_BINDING_ENERGYS[end] + VAC_MIGRATE_BARRIERS[1]
            end
            frequency = k0*exp(-barrier/universe.temperature/k_B)
            push!(universe.constants.vacEmitFrequencies, frequency)
        end
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
        frequencies = SiaFrequencies(universe, defect.size)
    else
        types = [1, 2]
        frequencies = VacFrequencies(universe, defect.size)
    end
    behaviors = Behaviors()
    behaviors.types = types
    behaviors.frequencies = frequencies
    frequency = sum(frequencies)
    behaviors.frequency = frequency
    defect.behaviors = behaviors
    push!(universe.defectFrequencies, frequency)
    universe.totalFrequency += frequency
end

function ChangeBehaviors!(universe::Universe, defect::Defect)
    oldFrequency = defect.behaviors.frequency
    if defect.type == 1
        frequencies = SiaFrequencies(universe, defect.size)
    else
        frequencies = VacFrequencies(universe, defect.size)
    end
    defect.behaviors.frequencies = frequencies
    frequency = sum(frequencies)
    defect.behaviors.frequency = frequency
    i = findfirst(x->x.index===defect.index, universe.defects)
    universe.defectFrequencies[i] = frequency
    universe.totalFrequency += frequency - oldFrequency
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
    deleteat!(universe.defectFrequencies, i)
    universe.totalFrequency -= defect.behaviors.frequency
end


function Base.delete!(cell::Cell, defect::Defect)
    deleteat!(cell.defects, findfirst(x->x.index==defect.index, cell.defects))
end

function Move!(universe::Universe, defect::Defect, coord::Vector{Float64})
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
    weights = Weights(universe.defectFrequencies)
    sample(universe.defects, weights)
end

function RandomABehavior(defect::Defect)
    weights = Weights(defect.behaviors.frequencies)
    sample(defect.behaviors.types, weights)
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
    r = rand()
    if defect.type == 1
        if r <= SIA_DISAPPEAR_RATE
            universe.record.sinkedSiaNum += defect.size
            delete!(universe, defect)
            return
        end
        sign = sample([1,-1])
        displace = SIA_DIRECTIONS[defect.directionIndex] * sign
        newCoord = defect.coord + displace
        Move!(universe, defect, newCoord)
    else
        if r <= VAC_DISAPPEAR_RATE
            universe.record.sinkedVacNum += defect.size
            delete!(universe, defect)
            return
        end
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
    ForwardTime!(universe)
    defect = RandomADefect(universe)
    type = RandomABehavior(defect)
    Behave!(universe, defect, type)
end

function ForwardTime!(universe::Universe)
    universe.time += -1/universe.totalFrequency*log(1-rand())
    # or simple verstion
    # universe.time += 1/universe.totalFrequency
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
    InitLog(logName)
    RefreshFile(dumpName)
    InitCells!(universe)
    InitRadius!(universe)
    universe.temperature = TEMPERATURE
    InitFrequencies!(universe)
    run(`tput sc`)
end


function Finish!(universe::Universe)
    run(`tput sc`)
end

function Test!(universe::Universe)
    RefreshFile(dumpName)
    InitLog(logName)
    Init!(universe)
    defect = Defect(rand(0:99,3), 2, rand(1:4), rand(10:20))
    push!(universe, defect)
    Dump(universe, dumpName)
end



function Run_small!(universe::Universe)
    Init!(universe)
    while universe.nStep <= 10000
        #exit()
        if universe.nStep & 10 == 0
            defect = Defect(rand(0.:199.,3), 1, rand(1:4), rand(1:1))
            push!(universe, defect)
        end
        IterStep!(universe)
        RecordSV!(universe)
        #print(universe)
        if universe.nStep & 100 == 0 
            println(universe.nStep) 
        end
        if universe.nStep & 100 == 0 
            Dump(universe, dumpName) 
        end
    end
end


# todo: 
# fix boundary cells ??????
# realistic frequency ??????
# beatifify screen output ??????
# outpot dataframe for python plot ??????
# include cascade data ???
