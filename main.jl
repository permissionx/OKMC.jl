using StatsBase
using Random
using Distributions

const SIA_DIRECTIONS = ([[1,1,1], 
                         [1,1,-1], 
                         [1,-1,1], 
                         [1,-1,-1]])

const FIRST_NEIGHBORS = [[1,1,1], 
                         [1,1,-1], 
                         [1,-1,1], 
                         [1,-1,-1],
                         [-1,1,1], 
                         [-1,1,-1], 
                         [-1,-1,1], 
                         [-1,-1,-1]]

const SECOND_NEIGHBORS = [[2,0,0], 
                          [-2,0,0], 
                          [0,2,0], 
                          [0,-2,0], 
                          [0,0,2], 
                          [0,0,-2]]

const SIA_DISAPPEAR_PROB = 1E-8


mutable struct Behaviors
    types::Vector{Float64}  #1 for diffusion, 2 for redirection
    probabilities::Vector{Float64}
    probability::Float64
    function Behaviors()
        types = Int64[]
        probabilities = Float64[]
        probability = 0.0
        new(types, probabilities, probability)
    end
end

mutable struct Defect
    index::Int64
    type::Int64
    coord::Vector{Int64}
    directionIndex::Int64
    size::Int64
    radius::Float64
    cellIndex::Int64
    behaviors::Behaviors
    function Defect(coord::Vector{Int64}, type::Int64, directionIndex::Int64, size::Int64)
        index = 0 
        cellIndex = 0
        behaviors = Behaviors()
        new(index, type, coord, directionIndex, size, 0., cellIndex, behaviors)
    end
end


mutable struct Cell
    index::Int64
    neighbors::Vector{Cell}
    defects::Vector{Defect}
    function Cell(index::Int64)
        new(index, Cell[], Defect[])
    end
end

mutable struct Universe
    nStep::Int64
    maxIndex::Int64
    mapSize::Vector{Int64}
    cellLength::Int64
    nsCells::Vector{Int64} #numbers of cells in 3 dimensions
    cells::Array{Cell, 3}
    defects::Vector{Defect}
    defectProbabilities::Vector{Float64}
    totalProbability::Float64
    function Universe(mapSize::Vector{Int64}, cellLength::Int64)
        nsCells = floor.(Int64, mapSize / cellLength)
        cells = Array{Cell, 3}(undef, nsCells[1], nsCells[2], nsCells[3])
        cellIndex = 0
        for i in 1:length(cells)
            cellIndex += 1
            cells[i] = Cell(cellIndex)
        end
        defects = Defect[]
        defectProbabilities = Float64[]
        maxIndex = 0
        totalProbability = 0 
        nStep = 0
        new(nStep, maxIndex, mapSize, cellLength, nsCells, cells, defects, defectProbabilities, totalProbability)
    end
end

CellCoord(universe::Universe, coord::Vector{Int64}) = floor.(Int64, coord/universe.cellLength .+ 1)

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
    write(file, "ITEM: ATOMS id type x y z size direction radius\n")
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


    

function GetCell(universe::Universe, cellCoord::Vector{Int64})
    for i in 1:3
        if cellCoord[i] > universe.nsCells[i]
            cellCoord[i] -= 1
        end
    end
    return universe.cells[cellCoord[1], cellCoord[2], cellCoord[3]]
end

function PBCCellCoord!(universe::Universe, cellCoord::Vector{Int64})
    for i in 1:3
        if cellCoord[i] < 1
            cellCoord[i] = universe.nsCells[i] + cellCoord[i]
        elseif cellCoord[i] > universe.nsCells[i]
            cellCoord[i] = cellCoord[i] - universe.nsCells[i]
        end
    end
end

function InitCells(universe::Universe)
    for i in 1:universe.nsCells[1]
        for j in 1:universe.nsCells[2]
            for k in 1:universe.nsCells[3]
                cell = universe.cells[i, j, k]
                for di in -1:1
                    for dj in -1:1
                        for dk in -1:1
                            neighborCellCoord = [i, j, k] + [di, dj, dk]
                            PBCCellCoord!(universe, neighborCellCoord)
                            neighborCell = GetCell(universe, neighborCellCoord)
                            push!(cell.neighbors, neighborCell)
                        end
                    end
                end
            end
        end
    end
end

function Base.display(cell::Cell)
    println("Cell: ", cell.index)
    println("id type coord directionIndex size radius cellIndex")
    for defect in cell.defects
        println("$(defect.index) $(defect.type) $(defect.coord) $(defect.directionIndex) \
        $(defect.size) $(defect.radius) $(defect.cellIndex)")
    end
    println("----------------------")
end

function Base.display(cells::Array{Cell})
    for cell in cells
        display(cell)
    end
end



function CrossCells!(universe::Universe, defect::Defect)
    cellCoord = CellCoord(universe, defect.coord)
    cell = GetCell(universe, cellCoord)
    if cell.index != defect.cellIndex
        delete!(universe.cells[defect.cellIndex], defect)
        push!(cell, defect)
    end
end

function Delta(universe::Universe, coord1::Vector{Int64}, coord2::Vector{Int64}) 
    # vactor points from 1 to 2 
    delta = Vector{Int64}(coord2 - coord1)
    for i in 1:3
        if delta[i] < -universe.mapSize[i] / 2
            delta[i] += universe.mapSize[i]
        elseif delta[i] > universe.mapSize[i] / 2
            delta[i] -= universe.mapSize[i]
        end
    end
    delta
end

function Distance(universe::Universe, coord1::Vector{Int64}, coord2::Vector{Int64})
    delta = Delta(universe, coord1, coord2)
    sqrt(sum(delta .* delta))
end

function FindNeighbors(universe::Universe, defect::Defect)
    cell = universe.cells[defect.cellIndex]
    neighbors = Defect[]
    for cell in cell.neighbors
        for neighbor in cell.defects
            if neighbor.index != defect.index
                if Distance(universe, neighbor.coord, defect.coord) <= neighbor.radius+defect.radius
                    push!(neighbors, neighbor)
                end
            end
        end
    end
    neighbors
end

function CoordInBCC(coord::Vector{Float64})
    # To do: should shift to the nearest lattice point, but does not matter seriously
    newCoord = Vector{Int64}(undef, 3)
    newCoord[3] = Int64(round(coord[3]))
    if iseven(newCoord[3])
        for i in 1:2
            x = floor(coord[i])
            if iseven(x)
                newCoord[i] = Int64(x)
            else
                newCoord[i] = Int64(x+1)
            end
        end
    else
        for i in 1:2
            x = floor(coord[i])
            if isodd(x)
                newCoord[i] = Int64(x)
            else
                newCoord[i] = Int64(x+1)
            end
        end
    end
    newCoord
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
        #newCoord = CoordInBCC(center)
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


function Move!(universe::Universe, defect::Defect, coord::Vector{Int64})
    PBCCoord!(universe, coord)
    defect.coord = coord
    CrossCells!(universe, defect)
    Changed!(universe, defect)
end


function Changed!(universe::Universe, defect::Defect)
    neighbors = FindNeighbors(universe, defect)
    if length(neighbors) > 0 
        neighbor = sample(neighbors)
        Reaction!(universe, defect, neighbor)
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


function SiaProbabilities(size::Int64)
    return [1.,0.001.]
end

function VacProbabilities(size::Int64)
    if size == 1
        return Float64[1.,0.5]
    else
        return Float64[0,0]
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


function Base.push!(universe::Universe, defect::Defect)
    cellCoord = CellCoord(universe, defect.coord)
    cell = GetCell(universe, cellCoord)
    push!(cell, defect)
    universe.maxIndex += 1
    defect.index = universe.maxIndex
    Radius!(defect)
    InitBehaviors!(universe, defect)
    push!(universe.defects, defect)
    Changed!(universe, defect)
end


function Base.push!(cell::Cell, defect::Defect)
    push!(cell.defects, defect)
    defect.cellIndex = cell.index
end

function Radius!(defect::Defect)
    defect.radius = ((3*defect.size)/(4*pi))^(1/3)
end

function ChangeSize!(universe::Universe, defect::Defect, size::Int64)
    defect.size = size
    Radius!(defect)
    ChangeBehaviors!(universe, defect)
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

function RandomADefect(universe::Universe)
    weights = Weights(universe.defectProbabilities)
    sample(universe.defects, weights)
end

function RandomABehavior(defect::Defect)
    sample(defect.behaviors.types, Weights(defect.behaviors.probabilities))
end

function Behave!(universe::Universe, defect::Defect, type::Int64)
    if type == 1
        MigrateFirstNeighbor!(universe, defect)
    elseif type == 2
        MigrateSecondNeighbor!(universe, defect)
    else
        Steer!(defect)
    end
end

function MigrateFirstNeighbor!(universe::Universe, defect::Defect)
    if defect.type == 1
        r = rand()
        if r <= SIA_DISAPPEAR_PROB
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
    Behave!(univesre, defect, type)
end

function Run!(universe::Universe)
    InitCells(universe)
    while universe.nStep < 1000000
        if universe.nStep % 1000 == 0
            defect = Defect(rand(0:99,3), rand(1:2), rand(1:4), rand(10:20))
            push!(universe, defect)
        end
        IterStep!(universe)
        if universe.nStep % 10000 == 0
            println(universe.nStep)
            Dump(universe, dumpName)
        end
    end
end

Random.seed!(31415926)
mapSize = [100,100,100]
cellLength = 20
universe = Universe(mapSize, cellLength)
const dumpName = "/mnt/c/Users/XUKE/Desktop/run.dump"
RefreshFile!(dumpName)
Run!(universe)

#empty
#using Profile, PProf
#@profile Run!(universe)
#pprof(;webport=58599)

