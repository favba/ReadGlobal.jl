using ReadGlobal
using Base.Test

for t in (Float32,Float64)
  for l in (1,2)
    a = rand(t,6,6,8)
    b = @view a[1:(end-l),:,:]
    write("test.test",a)

    @test b == readpadded("test.test",t,(6-l,6,8))
  end
end

rm("test.test")