using PyCall
@pyimport matplotlib.animation as animation
using PyPlot
using JLD


data = load("simulation_results_rungekutta.jld")
println("data loaded")
function save_movie(array, filename)
    fig = figure()
    ims = PyCall.PyObject([])
    for i in 1:10:size(array, 3)
        tmp_im = imshow(array[:,:,i], animated=true)
        push!(ims, PyCall.PyObject[tmp_im])
    end

    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=true, repeat_delay=1000)
    ani[:save](filename)
    close(fig)
end

save_movie(data["phi"], "phi.mp4")
save_movie(data["V"], "V.mp4")
save_movie(data["U"], "U.mp4")


println("finished")
