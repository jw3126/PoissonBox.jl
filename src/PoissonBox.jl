module PoissonBox

using TiledIteration
import LinearAlgebra
import IterativeSolvers

function CartesianIndices_inner(x::AbstractArray)
    ranges = map(axes(x)) do ax
        (first(ax)+1):(last(ax)-1)
    end
    CartesianIndices(ranges)
end

function forward!(out,u)
    if axes(out) != axes(u)
        msg = """
        Axes do not match:
        axes(out) == $(axes(out))
        axes(u) == $(axes(u))
        """
        throw(DimensionMismatch(msg))
    end
    all_inds = CartesianIndices(u)
    inner_inds = CartesianIndices_inner(u)
    begin
        for I in inner_inds
            out[I] = laplace(u, I)
        end
        for I in EdgeIterator(all_inds, inner_inds)
            out[I] = u[I]
        end
    end
    return out
end
function forward(u)
    out = similar(u, float(eltype(u)))
    forward!(out, u)
end

Base.@propagate_inbounds function laplace(u::AbstractArray{<:Any,1},I::CartesianIndex)
    (i,) = Tuple(I)
    (u[i+1] + u[i-1])/2 - u[i]
end

Base.@propagate_inbounds function laplace(u::AbstractArray{<:Any,2},I::CartesianIndex)
    (i,j) = Tuple(I)
    (u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i, j-1])/4 - u[i,j]
end
Base.@propagate_inbounds function laplace(u::AbstractArray{<:Any,3},I::CartesianIndex)
    (i,j,k) = Tuple(I)
    (u[i+1,j  ,k  ] + u[i-1,j   ,k  ] + 
     u[i  ,j+1,k  ] + u[i  , j-1,k  ] +
     u[i  ,j  ,k+1] + u[i  , j  ,k-1]
    )/6 - u[i,j,k]
end

################################################################################
#### LaplaceDirichlet
################################################################################
struct LaplaceDirichlet{T,N}
    size::NTuple{N,Int}
end

make_eltype_float(x::AbstractArray{<:AbstractFloat}) = x
function make_eltype_float(x::AbstractArray)
    T = float(eltype(x))
    out = similar(x, T)
    copy!(out, x)
end


function poisson_problem(arr::AbstractArray{T,N}) where {T,N}
    op = LaplaceDirichlet{float(T),N}(size(arr))
    v = reshape(make_eltype_float(arr), length(arr))
    op, v
end
function reshape_solution(op, sol)
    reshape(sol, op.size)
end

function LinearAlgebra.mul!(out, op::LaplaceDirichlet, v::AbstractVector)
    forward!(reshape(out, op.size), reshape(v, op.size))
    out
end
function Base.:(*)(op::LaplaceDirichlet, v::AbstractVector)
    T = typeof(one(eltype(v))/one(eltype(op)))
    out = similar(v,T)
    mul!(out, op, v)
end
function Base.eltype(::Type{LaplaceDirichlet{T,N}}) where {T,N}
    T
end
function Base.size(op::LaplaceDirichlet, d)
    (1,2)[d] # for appropriate error
    prod(op.size)
end

function Base.collect(::Type{T}, op::LaplaceDirichlet) where {T}
    # useful for debugging
    n = size(op,1)
    out = Matrix{T}(undef, n,n)
    v = zeros(Int,n)
    for i in 1:n
        v[i] = 1
        if i > 1
            v[i-1] = 0
        end
        mul!(view(out,:,i), op, v)
    end
    out
end
function Base.collect(op::LaplaceDirichlet)
    Base.collect(eltype(op), op)
end

function solve(solver, y::AbstractArray, solver_args...; solver_kw...)
    op, v = poisson_problem(y)
    sol = solver(op, v, solver_args...; solver_kw...)
    reshape_solution(op, sol)
end
function solve(y::AbstractArray)
    solve(IterativeSolvers.bicgstabl, y, 2)
end

end

