__precompile__()
module ReadGlobal

using DelimitedFiles

export readglobal, getdimsize, readpadded, readpadded!, readfield, readfield!, checkinput, getnfilter, doinchunks, read_info

function findglobal()
    filename="global"
    if !isfile(filename)
        for i = 1:10
            filename = "../" * filename
            isfile(filename) && break
        end
    end
    isfile(filename) || throw(error("Could not find global file"))
    return readdlm(filename,String,comment_char='/')
end

function readglobal()
    text = findglobal()
    dict = Dict{Symbol,String}()
    for i in 1:size(text)[1]
        dict[Symbol(text[i,1])] = text[i,2]
    end
    return dict
end

function getdimsize()
    dict = readglobal()
    nx = parse(Int,dict[:nx])
    ny = parse(Int,dict[:ny])
    nz = parse(Int,dict[:nz])
    x  = parse(Float64,dict[:xDomainSize])
    y  = parse(Float64,dict[:yDomainSize])
    z  = parse(Float64,dict[:zDomainSize])

    return nx,ny,nz,x,y,z
end

function readpadded!(stream::IO, field::AbstractArray{T,3}) where {T}
    dim = size(field)
    nx,ny,nz = dim
    sb = sizeof(T)
    npad = iseven(nx) ? 2 : 1
    ind = LinearIndices(dim)
    for k in 1:nz
        for j in 1:ny
            unsafe_read(stream,Ref(field,ind[1,j,k]),nx*sb)
            skip(stream, npad*sb)
        end
    end
    return field
end

function readpadded!(file::AbstractString, field::AbstractArray{T,N}) where {T,N}
    open(file) do io
        dims = size(field)
        ndim = [dims...]
        ndim[1] += iseven(ndim[1]) ? 2 : 1
        filesize(file) !== prod(ndim)*sizeof(T) && error("Dimension Mismatch")
        return readpadded!(io,field) 
    end
end

function readpadded(stream,T::DataType,dims)
    field = Array{T}(undef,dims)
    return readpadded!(stream,field)
end

function readpadded(stream,dims)
    field = Array{Float64}(undef,dims)
    return readpadded!(stream,field)
end

function checkinput(filename::String,nx::Int,ny::Int,nz::Int)
    sizefile = filesize(filename)
    if sizefile == (nx+2)*ny*nz*8
        dtype = Float64
        padded = true
    elseif sizefile == nx*ny*nz*8
        dtype = Float64
        padded = false
    elseif sizefile == (nx+2)*ny*nz*4
        dtype = Float32
        padded = true
    elseif sizefile == nx*ny*nz*4
        dtype = Float32
        padded = false
    else
        error("""Input file "$(filename)" not matching number of gridpoints on global file (nx:$nx,ny:$ny,nz:$nz)""")
    end
    return dtype,padded
end

function readfield(filename::String,nx::Int,ny::Int,nz::Int)
    dtype,padded = checkinput(filename,nx,ny,nz)
    return padded ? readpadded(filename,dtype,(nx,ny,nz)) : read!(filename,Array{dtype}(undef,nx,ny,nz))
end

function readfield(filename::String)
    nx,ny,nz,x,y,z = getdimsize()
    return readfield(filename,nx,ny,nz)
end

function readfield!(filename::String,field::DenseArray)
    dtype,padded = checkinput(filename,size(field)...)
    return padded ? readpadded!(filename,field) : read!(filename,field)
end

function getnfilter()
    place = split(pwd(),"/")
    if occursin(r"^Model_",place[end])
        N = place[end-1][2:end]
        Fil = place[end-2][1]
    else
        N = place[end][2:end]
        Fil = place[end-1][1]
    end
    return N, Fil
end

function best_nc(ni::Int,no::Int,nx::Int,ny::Int,nz::Int)
    nc::Int = 1
    totalmem::Int = floor(Sys.total_memory())
    use::Float64 = 0.75
    usemem = use*totalmem
    totalfields::Int = (ni+no)*nx*ny*nz*sizeof(Float64)
    if usemem < totalfields
        newsize = totalfields
        while usemem < newsize
            nc *= 2
            newsize /= 2
        end
    end
    return nc
end

init(T,L) = T(L...)
init(::Type{Array{T,N}},L) where {T,N} = Array{T,N}(undef,L...)

function doinchunks(func::Function,nc::Int=0;input::NTuple{Ni,Pair}=(),output::NTuple{No,Pair}=()) where {Ni,No}

    nx, ny, nz, xs, ys, zs = getdimsize()

    TNi = sum(length.(getfield.(input,:second)))
    TNo = sum(length.(getfield.(output,:second)))

    nc == 0 && (nc = best_nc(TNi,TNo,nx,ny,nz))

    ifnames = getfield.(input,:second)
    ofnames = getfield.(output,:second)

    @assert mod(nx*ny*nz,nc)==0
    chunk = Int(nx*ny*nz/nc)

    iarrays = init.(getfield.(input,:first),chunk) 
    oarrays = init.(getfield.(output,:first),chunk) 

    calculation(func,nc,iarrays,oarrays,ifnames,ofnames)

    return 0
end

mywrite(f,v) = write(f,v)
myread!(f,v) = read!(f,v)
mywrite(f::Tuple{<:Union{<:IO,<:AbstractString}},v::Array{T,N}) where {T,N} = write(f[1],v)
myread!(f::Tuple{<:Union{<:IO,<:AbstractString}},v::Array{T,N}) where {T,N} = read!(f[1],v)

function calculation(func::F,nc,iarrays,oarrays,ifn,ofn) where {F<:Function}

    ifiles = (x->open.(x,"r")).(ifn)
    ofiles = (x->open.(x,"w")).(ofn)

    for j=1:nc
        for i in Base.OneTo(length(iarrays))
            myread!(ifiles[i],iarrays[i])
        end

        func(iarrays...,oarrays...)

        for i in Base.OneTo(length(oarrays))
            mywrite(ofiles[i],oarrays[i])
        end

        #gc()

    end

    for file in ifiles
        close.(file)
    end

    for file in ofiles
        close.(file)
    end

end

function read_info(stream)
    med = parse(Float64,readline(stream)[7:end])
    std = parse(Float64,readline(stream)[5:end])
    maxv = parse(Float64,split(readline(stream))[2])
    minv = parse(Float64,split(readline(stream))[2])
    return med,std,maxv,minv
end

read_info(stream::AbstractString) = (open(stream) do f; return read_info(f);end)

end # module
