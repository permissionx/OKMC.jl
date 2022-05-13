using StatsBase

abstract type Defect end

mutable struct Behaviors
    displaces::Vector{Vector{Int64}}
    probabilities::Vector{Float64}
end

mutable struct Sia <: Defect
    index::Int64
    coord::Vector{Int64}
    directionIndex::Int64
    size::Int64
    radius::Float64
    cellIndex::Int64
    function Sia(coord::Vector{Int64}, directionIndex::Int64, size::Int64)
        index = 0 
        cellIndex = 0
        new(index, coord, directionIndex, size, 0., cellIndex)
    end
end

mutable struct Vac <: Defect
    index::Int64
    coord::Vector{Int64}
    directionIndex::Int64
    size::Int64
    radius::Float64
    cellIndex::Int64
    function Vac(coord::Vector{Int64}, directionIndex::Int64, size::Int64)
        index = 0
        cellIndex = 0
        new(index, coord, directionIndex, size, 0., cellIndex)
    end
end


mutable struct Cell
    index::Int64
    neighbors::Vector{Cell}
    sias::Vector{Sia}
    vacs::Vector{Vac}
    function Cell(index::Int64)
        new(index, Cell[], Sia[], Vac[])
    end
end

mutable struct Universe
    maxIndex::Int64
    mapSize::Vector{Int64}
    cellLength::Int64
    nsCells::Vector{Int64} #numbers of cells in 3 dimensions
    cells::Array{Cell, 3}
    sias::Vector{Sia}
    vacs::Vector{Vac}
    function Universe(mapSize::Vector{Int64}, cellLength::Int64)
        nsCells = floor.(Int64, mapSize / cellLength)
        cells = Array{Cell, 3}(undef, nsCells[1], nsCells[2], nsCells[3])
        cellIndex = 0
        for i in 1:length(cells)
            cellIndex += 1
            cells[i] = Cell(cellIndex)
        end
        new(0, mapSize, cellLength, nsCells, cells, Sia, Vac[])
    end
end

function NDefects(universe::Universe)
    nSia = 0
    nVac = 0
    for cell in universe.cells
        nSia += length(cell.sias)
        nVac += length(cell.vacs)
    end
    nSia, nVac, nSia + nVac
end

CellCoord(universe::Universe, coord::Vector{Int64}) = floor.(Int64, coord/universe.cellLength .+ 1)


function GetCell(universe::Universe, cellCoord::Vector{Int64})
    try
        return universe.cells[cellCoord[1], cellCoord[2], cellCoord[3]]
    catch BoundsError
        return universe.cells[cellCoord[1]-1, cellCoord[2]-1, cellCoord[3]-1]
    end
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
    print("Sias: ")
    for sia in cell.sias
        print("$(sia.index) ")
    end
    print("\n")
    print("Vacs: ")
    for vac in cell.vacs
        print("$(vac.index) ")
    end
    print("\n----------------------\n")
end

function Base.display(cells::Array{Cell})
    for cell in cells
        display(cell)
    end
end

function Base.push!(universe::Universe, defect::Defect)
    cellCoord = CellCoord(universe, defect.coord)
    PBCCellCoord!(universe, cellCoord)
    cell = GetCell(universe, cellCoord)
    push!(cell, defect)
    universe.maxIndex += 1
    defect.index = universe.maxIndex
    Radius!(defect)
    Changed!(universe, defect)
end

function CrossCells(universe::Universe, defect::Defect)
    cellCoord = CellCoord(universe, defect.coord)
    PBCCellCoord!(universe, cellCoord)
    cell = GetCell(universe, cellCoord)
    if cell.index != defect.cellIndex
        delete!(cell, defect)
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
    cell = defect.cell
    neighborSias = Sia[]
    neighborVacs = Vac[]
    for neighbor in cell.neighbors
        for sia in neighbor.sias
            if sia.directionIndex != defect.directionIndex
                if Distance(universe, sia.coord, defect.coord) <= sia.radius+defect.radius
                    push!(neighborSias, sia)
                end
            end
        end
        for vac in neighbor.vacs
            if vac.directionIndex != defect.directionIndex
                if Distance(universe, vac.coord, defect.coord) <= vac.radius+defect.radius
                    push!(neighborVacs, vac)
                end
            end
        end
    end
    neighborSias, neighborVacs
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

function Reaction!(universe::Universe, defect1::T, defect2::T) where T<:Union{Vac, Sia}
    center = (defect1.coord*defect.size + defect2.coord*defect.size) / (defect1.size + defec2.size)
    newCoord = CoordInBCC(center)
    PBCCellCoord!(universe, newCoord)
    defect = (defect1.size > defect2.size) ? defect1 : defect2
    setSize!(defect, defect1.size + defect2.size)
    delete!(universe, defect2)
    Move!(defect, newCoord)
end

function Reaction!(universe::Universe, defect1::Union{Vac, Sia}, defect2::Union{Vac, Sia})
    if defect1.size == defect2.size
        delete!(universe, defect1)
        delete!(universe, defect2)
        return
    elseif defect1.size > defect2.size
        Swallow!(universe, defect1, defect2)
    else
        Swallow!(universe, defect2, defect1)
    end
end

function Swallow!(universe::Universe, defect1::Union{Vac, Sia}, defect2::Union{Vac, Sia})
    radius = defect1.radius 
    SetSize!(defect1, defect1 - defect2)
    moveLength = radius - defect1.radius
    delta = Delta(universe, defect1.coord, defect2.coord)
    distance = sqrt(sum(delta .* delta))
    newCoordFloat = Vector{Float64}(undef,3)
    for i in 1:3
        newCoordFloat[i] -= moveLength * delte[i] / distance
    end
    newCoord = CoordInBCC(newCoordFloat)
    PBCCellCoord!(universe, newCoord)
    delete!(universe, defect2)
    Move!(defect1, newCoord)
end
    

function Move!(defect::Defect, coord::Vector{Int64})
    defect.coord = coord
    CrossCells(universe, defect)
    Changed!(universe, defect)
end


function Changed!(universe::Universe, defect::Defect)
    sias, vacs = FindNeighbors(universe, defect)
    if length(sias) == 0 && length(vacs) == 0
        return
    elseif length(vacs) == 0
        sia = sample(sias)
        Reaction!(universe, defect, sia)
    elseif length(sias) == 0
        vac = sample(vacs)
        Reaction!(universe, defect, vac)
    else
        typeNum = rand(1:length(sias)+length(vacs))
        if typeNum <= length(sias)
            sia = sample(sias)
            Reaction!(universe, defect, sia)
        else
            vac = sample(vacs)
            Reaction!(universe, defect, vac)
        end
    end
end


function Base.delete!(cell::Cell, sia::Sia)
    deleteat(cell.sias, findfirst(x->x.index==sia.index, cell.sias))
end

function Base.delete!(cell::Cell, vac::Vac)
    deleteat(cell.vacs, findfirst(x->x.index==vac.index, cell.vacs))
end

function Base.delete!(universe::Universe, defect::Defect)
    cellCoord = CellCoord(universe, defect.coord)
    PBCCellCoord!(universe, cellCoord)
    cell = GetCell(universe, cellCoord)
    delete!(cell, defect)
end

function Base.push!(cell::Cell, sia::Sia)
    push!(cell.sias, sia)
    sia.cellIndex = cell.index
end

function Base.push!(cell::Cell, vac::Vac)
    push!(cell.vacs, vac)
    vac.cellIndex = cell.index
end

function Radius!(defect::Defect)
    defect.radius = ((3*defect.size)/(4*pi))^(1/3)
end

function SetSize!(defect::Defect, size::Int64)
    defect.size = size
    Radius!(defect)
end



mapSize = [100,100,100]
cellLength = 20
universe = Universe(mapSize, cellLength)
InitCells(universe)
sia = Sia([50,50,50], 1, 10)
vac = Vac([40,40,40], 1, 12)
push!(universe, sia)
push!(universe, vac)
