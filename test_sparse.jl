function trace_sparse1!(A  :: SparseMatrixCSC{Float64,Int64}, ind :: Array{Int,2}, row :: Array{Int,1}, col :: Array{Int,1}, lin :: Array{Int,1})

    for k in 1:size(ind,1)
        A[ind[k,:],ind[k,:]] *= 10.0
    end

    return nothing
end


function trace_sparse2!(A  :: SparseMatrixCSC{Float64,Int64}, ind :: Array{Int,2}, row :: Array{Int,1}, col :: Array{Int,1}, lin :: Array{Int,1})

    for k in 1:size(row,1)
        A[row[k],col[k]] *= 2.0
    end

    return nothing
end


function trace_sparse3!(A  :: SparseMatrixCSC{Float64,Int64}, ind :: Array{Int,2}, row :: Array{Int,1}, col :: Array{Int,1}, lin :: Array{Int,1})

    for k in lin
        A[k] *= 2.0
    end

    return nothing
end

function trace_sparse4!(A  :: SparseMatrixCSC{Float64,Int64}, ind :: Array{Int,2}, row :: Array{Int,1}, col :: Array{Int,1}, lin :: Array{Int,1})

    for k in 1:nnz(A)
        A[k] *= 2.0
    end

    return nothing
end

# -------------------------------------------------------
# -------------------------------------------------------
# -------------------------------------------------------

function create_ind(n :: Int)

    ind = rand(1:n, n, 2)
    row = vec(ind[:,[1;1;2;2]]')
    col = vec([ind ind]')
    
    return ind, row, col
end

function create_array(ind :: Array{Int,2}, row :: Array{Int,1}, col :: Array{Int,1})

    A = sparse(row, col, zeros(Float64, length(row)), length(row), length(row))
    
    return A
end

# -------------------------------------------------------
# -------------------------------------------------------
# -------------------------------------------------------

function fill_sparray(ind :: Array{Int,2}, row :: Array{Int,1}, col :: Array{Int,1})

    val = fill_val(size(ind,1))

    A = sparse(row, col, val[:], length(row), length(row))

    return A
end

function fill_val(n :: Int)

    val = Array{Float64,3}(2,2,n)
    
    for i=1:n
        val[:,:,i] = rand(2,2)
    end

    return val
end

# -------------------------------------------------------
# -------------------------------------------------------
# -------------------------------------------------------

function fill_sparray!(A :: SparseMatrixCSC{Float64,Int64}, ind :: Array{Int,2}, row :: Array{Int,1}, col :: Array{Int,1})

    sub = sub2ind(size(A), row, col)
    
    for i in sub
        A[i] += rand()
        #fill_val!(view(A,ind[i,:], ind[i,:]))
    end

    return nothing
end

function fill_val!(val_loc :: SubArray{Float64,2,SparseMatrixCSC{Float64,Int64},Tuple{Array{Int64,1},Array{Int64,1}},false})
    
        @time val_loc += rand(2,2)

    return nothing
end
