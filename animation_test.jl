using PyCall
@pyimport matplotlib.animation as animation
using PyPlot


A = randn(20,20,20)

fig = figure()
ims = PyCall.PyObject([])
for i in 1:20
    tmp_im = imshow(A[:,:,i], animated=true)
    push!(ims, PyCall.PyObject[tmp_im])
end

ani = animation.ArtistAnimation(fig, ims, interval=50, blit=true, repeat_delay=1000)

ani[:save]("test.mp4")

println("finished")
