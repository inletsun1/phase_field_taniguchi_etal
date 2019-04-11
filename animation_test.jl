using PyCall
matplotlib = pyimport("matplotlib")
matplotlib.use("tkagg")
animation = pyimport("matplotlib.animation")
plt = pyimport("matplotlib.pyplot")


A = randn(20,20,20)

fig = plt.figure()
ims = PyCall.PyObject([])
for i in 1:20
    tmp_im = plt.imshow(A[:,:,i], animated=true)
    push!(ims, PyCall.PyObject[tmp_im])
end

println("hoge")
ani = animation.ArtistAnimation(fig, ims, interval=50, blit=true, repeat_delay=1000)

ani.save("test.mp4", writer="html")

println("finished")
