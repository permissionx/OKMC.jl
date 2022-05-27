using StatsBase
using Random
using Distributions
using Crayons

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

mutable struct Constants
    siaRadii::Vector{Float64}
    vacRadii::Vector{Float64}
    vacMigrateFrequencies::Vector{Float64}
    siaMigrateFrequencies::Vector{Float64}
    siaSteerFrequency::Float64
    vacEmitFrequencies::Vector{Float64}
    function Constants()
        siaRadii = Float64[]
        vacRadii = Float64[]
        vacMigrateFrequencies = Float64[]
        siaMigrateFrequencies = Float64[]
        siaSteerFrequency = 0.0
        vacEmitFrequencies = Float64[]
        new(siaRadii, vacRadii, vacMigrateFrequencies, siaMigrateFrequencies, siaSteerFrequency, vacEmitFrequencies)
    end
end 

mutable struct Behaviors
    types::Vector{Int64}  #1 for diffusion, 2 for redirection
    frequencies::Vector{Float64}
    frequency::Float64
    function Behaviors()
        types = Int64[]
        frequencies = Float64[]
        frequency = 0.0
        new(types, frequencies, frequency)
    end
end

mutable struct Defect
    index::Int64
    type::Int64
    coord::Vector{Float64}
    directionIndex::Int64
    size::Int64
    radius::Float64
    cellIndex::Int64
    behaviors::Behaviors
    function Defect(coord::Vector{Float64}, type::Int64, directionIndex::Int64, size::Int64)
        index = 0 
        cellIndex = 0
        behaviors = Behaviors()
        new(index, type, coord, directionIndex, size, 0., cellIndex, behaviors)
    end
end

mutable struct Cell
    index::Int64
    normalNeighbors::Vector{Cell}
    crossNeighbors::Vector{Cell}
    crossSigns::Vector{Vector{Int64}}
    defects::Vector{Defect}
    function Cell(index::Int64)
        crossSigns = Vector{Vector{Int64}}()
        new(index, Cell[], Cell[], crossSigns, Defect[])
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


mutable struct History
    steps::Vector{Int64}
    times::Vector{Float64}
    siaNums::Vector{Int64}
    vacNums::Vector{Int64}
    siaClusterNums::Vector{Int64}
    vacClusterNums::Vector{Int64}
    sinkedSiaNums::Vector{Int64}
    annihilatedSiaNums::Vector{Int64}
    function History()
        steps = Int64[]
        times = Float64[]
        siaNums = Int64[]
        vacNums = Int64[]
        siaClusterNums = Int64[]
        vacClusterNums = Int64[]
        sinkedSiaNums = Int64[]
        annihilatedSiaNums = Int64[]
        new(steps, times, siaNums, vacNums, siaClusterNums, vacClusterNums, sinkedSiaNums, annihilatedSiaNums)
    end
end


mutable struct Record
    sinkedSiaNum::Int64
    annihilatedSiaNum::Int64
    function Record()
        new(0,0)
    end
end


mutable struct Universe
    nStep::Int64
    time::Float64
    maxIndex::Int64
    mapSize::Vector{Float64}
    cellLength::Int64
    nsCells::Vector{Int64} #numbers of cells in 3 dimensions
    cells::Array{Cell, 3}
    defects::Vector{Defect}
    defectFrequencies::Vector{Float64}
    totalFrequency::Float64
    constants::Constants
    history::History
    record::Record
    function Universe(mapSize::Vector{Float64}, cellLength::Int64)
        nsCells = floor.(Int64, mapSize / cellLength)
        cells = Array{Cell, 3}(undef, nsCells[1], nsCells[2], nsCells[3])
        cellIndex = 0
        for i in 1:length(cells)
            cellIndex += 1
            cells[i] = Cell(cellIndex)
        end
        defects = Defect[]
        defectFrequencies = Float64[]
        maxIndex = 0
        totalFrequency = 0 
        nStep = 0
        constants = Constants()
        history = History()
        record = Record()
        time = 0.0
        new(time, nStep, maxIndex, mapSize, cellLength, nsCells, cells, defects, 
            defectFrequencies, totalFrequency, constants, history, record)
    end
end
