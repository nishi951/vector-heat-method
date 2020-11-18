using LinearAlgebra


# Generate a flat triangle mesh with simple connectivity
function make_flatmesh(height, width, dx)
    N = height*width
    dy = dx * √(3)/2
    square_grid = [[dy*y, dx*x, 0] for y in 1:height, x in 1:width]

    offset = dx/2

    for row in 1:height
        if (row % 2) == 0
            for col in 1:width
                square_grid[row, col][2] += dx/2
            end
        end
    end

    # build faces
    F = Array{Array{Int64, 1},1}()
    lin = LinearIndices(square_grid)
    for row in 1:height-1
        if row % 2 == 1
            for col in 1:width-1
                push!(F, [lin[row, col], lin[row, col+1], lin[row+1, col]])
                push!(F, [lin[row+1, col], lin[row, col+1], lin[row+1, col+1]])
            end
        else
            for col in 1:width-1
                push!(F, [lin[row, col], lin[row, col+1], lin[row+1, col+1]])
                push!(F, [lin[row+1, col], lin[row, col], lin[row+1, col+1]])
            end
        end
    end

    V = Array(cat(square_grid[:]..., dims=2)')
    return V, F
end
