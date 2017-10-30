__precompile__()
module ReadGlobal

export readglobal, getdimsize, readpadded, readpadded!, readfield, checkinput, getnfilter, doinchunks

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

function readpadded!(stream::IO, field::AbstractArray{T,N}) where {T,N}
  if N == 1
    read!(stream,field)
  else
    dims = size(field)
    nx = dims[1]
    nb = sizeof(T)*nx
    npencils = prod(dims)/nx
    npad = iseven(nx)? 2 : 1
    for i=0:(npencils-1)
      unsafe_read(stream,Ref(field,Int((nx+npad)*i+1)),nb)
    end
  end
  return field
end

function readpadded!(file::AbstractString, field::AbstractArray)
  open(file) do io 
   return readpadded!(io,field) 
  end
end

function readpadded(stream,T::DataType,dims)
  field = Array{T}(dims)
  return readpadded!(stream,field)
end

function readpadded(stream,dims)
  field = Array{Float64}(dims)
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
  return padded ? readpadded(filename,dtype,(nx,ny,nz)) : read(filename,dtype,(nx,ny,nz))
end

function readfield(filename::String)
  nx,ny,nz,x,y,z = getdimsize()
  dtype,padded = checkinput(filename,nx,ny,nz)
  return padded ? readpadded(filename,dtype,(nx,ny,nz)) : read(filename,dtype,(nx,ny,nz))
end

function getnfilter()
  place = split(pwd(),"/")
  if ismatch(r"^Model_",place[end])
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
  use::Float64 = 0.65
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

function doinchunks(func::Function,nc::Int=0;input::NTuple{Ni,String}=(),output::NTuple{No,String}=()) where {Ni,No}

  nx, ny, nz, xs, ys, zs = getdimsize()

  nc == 0 && (nc = best_nc(Ni,No,nx,ny,nz))

  ifiles = [open(i,"r") for i in input]
  ofiles = [open(i,"w") for i in output]

  assert(mod(nx*ny*nz,nc)==0)
  chunk = Int(nx*ny*nz/nc)

  iarrays = [Array{Float64}((chunk,)) for i=1:Ni]
  oarrays = [Array{Float64}((chunk,)) for i=1:No]

  for j=1:nc

    for i in 1:Ni
      read!(ifiles[i],iarrays[i])
    end

    func(iarrays...,oarrays...)

    for i in 1:No
      write(ofiles[i],oarrays[i])
    end

    gc()

  end

  for file in ifiles
    close(file)
  end

  for file in ofiles
    close(file)
  end

return 0
end

end # module
