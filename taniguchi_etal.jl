using PyCall
@pyimport matplotlib.animation as animation
using PyPlot

#parameters
XSIZE = 100 + 2
YSIZE = 100 + 2
TSIZE = 1000 
dx = 0.5
dy = dx
Du = 0.05
Dv = 0.20
dt = 0.01

S = 1.0
tau = 0.83
eps = 1.0
eta = 1.0
M = 0.5
A0 = 5.0^2*pi
R = 5.0
#Fig. 4F
alpha = 8.0
beta = 4.8
Kk = 3.5
Kp = 5.2
gamma = 0.15
mu = 1.0
Xu = 50.0
Xv = 50.0

a = 10.0
b = 3.0e-2

println("dt = ", dt, ", 1/(2Du(1/dx^2 + 1/dy^2)) = ", 1/(2*Du*(1/(dx^2) + 1/(dy^2))))


#initialization
phi = zeros(Float64, YSIZE, XSIZE, TSIZE);
U = zeros(Float64, YSIZE, XSIZE, TSIZE);
V = zeros(Float64, YSIZE, XSIZE, TSIZE);
for i in 1:YSIZE
    for j in 1:XSIZE
        if sqrt((i-YSIZE/2)^2 + (j-XSIZE/2)^2) <= R/dx
            phi[i, j, 1] = 1.0
            V[i,j,1] = 1.0
            U[i,j,1] = 6.0

        end
    end
end

#boundary condition and diffusion discritization matrix
Ax = zeros(Float64, XSIZE, XSIZE)
Bx = zeros(Float64, XSIZE, XSIZE)
Ay = zeros(Float64, YSIZE, YSIZE)
By = zeros(Float64, YSIZE, YSIZE)
for i in 2:YSIZE-1
    Ay[i, i-1:i+1] = [1.0 -2.0 1.0] 
end
for i in 2:XSIZE-1
    Ax[i-1:i+1, i] = [1.0; -2.0; 1.0] 
end
for i in 2:YSIZE
    By[i, i-1:i] = [-1.0 1.0] 
end
for i in 2:XSIZE
    Bx[i-1:i, i] = [-1.0; 1.0] 
end
#=
#Dirichlet

Ay[1, 1:3] = [1.0 -2.0 1.0]
Ay[end, end-2:end] = [1.0 -2.0 1.0]
Ax[1:3, 1] = [1.0; -2.0; 1.0]
Ax[end-2:end, end] = [1.0; -2.0; 1.0]
=#


#Periodic
Ay[1, 1:2] = [-2.0 1.0]; Ay[1, end] = 1.0;
Ay[end, end-1:end] = [1.0 -2.0]; Ay[end, 1] = 1.0;
Ax[1:2, 1] = [-2.0; 1.0]; Ax[end, 1] = 1.0;
Ax[end-1:end, end] = [-2.0; 1.0]; Ax[1, end] = 1.0;
By[1, 1] = 1.0; By[1, end] = -1.0
Bx[1, 1] = 1.0; Bx[end, 1] = -1.0

println("start calculation")
#time development
for t in 2:TSIZE
    A = sum(phi[:,:,t-1].*dx.*dy)

    abs_nabra_phi = sqrt.(((1.0 ./dx) .* phi[:,:,t-1] *Bx).^2 .+ ((1.0 ./dy) .* By *phi[:,:,t-1]).^2)

    del_phi = (eta.*((1.0 ./dy) .^2.0 .* Ay * phi[:,:,t-1] .+ (1.0 ./dx) .^2.0 .* phi[:,:,t-1] * Ax .- 18.0.*phi[:,:,t-1].^2.0 .*(1.0 .- phi[:,:,t-1]).^2.0 ) .-M.*(A .- A0).*abs_nabra_phi .+ (a.*V[:,:,t-1] .- b.*U[:,:,t-1]).*abs_nabra_phi)./tau

    phi[:,:,t] = phi[:,:,t-1] .+ dt.*del_phi

    reaction_common = -alpha.*U[:,:,t-1].*V[:,:,t-1].^2 ./(Kk + sum(phi[:,:,t-1].*V[:,:,t-1].^2 .*dx.*dy)./A) .+beta.*U[:,:,t-1].*V[:,:,t-1]./(Kp.+sum(phi[:,:,t-1].*U[:,:,t-1].*dx.*dy)./A) 

    reaction_U = phi[:,:,t-1].*(reaction_common .+ S .-gamma.*U[:,:,t-1]).-Xu.*U[:,:,t-1].*abs_nabra_phi.^2 ./(sum(dx.*dy.*abs_nabra_phi.^2))

    reaction_V = phi[:,:,t-1].*(-reaction_common .-mu.*V[:,:,t-1]).-Xv.*V[:,:,t-1].*abs_nabra_phi.^2 ./(sum(dx.*dy.*abs_nabra_phi.^2))

    U[:,:,t] = U[:,:,t-1].+ dt .*(Du.*(phi[:,:,t-1] .*((1.0 ./dy) .^2 .* Ay * U[:,:,t-1] .+ (1.0 ./dx) .^2 .* U[:,:,t-1] * Ax) .+ ((1.0 ./dx) .* phi[:,:,t-1] *Bx .+ (1.0 ./dy) .* By *phi[:,:,t-1]) .*((1.0 ./dx) .* U[:,:,t-1] *Bx .+ (1.0 ./dy) .* By *U[:,:,t-1])) .+ reaction_U .- del_phi)

    V[:,:,t] = V[:,:,t-1].+ dt .*(Dv.*(phi[:,:,t-1] .*((1.0 ./dy) .^2 .* Ay * V[:,:,t-1] .+ (1.0 ./dx) .^2 .* V[:,:,t-1] * Ax) .+ ((1.0 ./dx) .* phi[:,:,t-1] *Bx .+ (1.0 ./dy) .* By *phi[:,:,t-1]) .*((1.0 ./dx) .* V[:,:,t-1] *Bx .+ (1.0 ./dy) .* By *V[:,:,t-1])) .+ reaction_V .- del_phi)
end
println("end calculation")

#make the movies
ims_phi = []
ims_U = []
ims_V = []

fig = figure()
for i=1:TSIZE
    #tmp_im = imshow(phi[:,:,i], animated=true)
    #push!(ims_phi, PyCall.PyObject[tmp_im])
    #tmp_im = imshow(U[:,:,i], animated=true)
    #push!(ims_U, PyCall.PyObject[tmp_im])
    tmp_im = imshow(V[:,:,i], animated=true)
    push!(ims_V, PyCall.PyObject[tmp_im])
end


#ani = animation.ArtistAnimation(fig, ims_phi, interval=100, blit=true, repeat_delay=1000)
#ani[:save]("phi.mp4")

#ani = animation.ArtistAnimation(fig, ims_U, interval=100, blit=true, repeat_delay=1000)
#ani[:save]("U.mp4")

ani = animation.ArtistAnimation(fig, ims_V, interval=100, blit=true, repeat_delay=1000)
ani[:save]("V.mp4")

println("finish")
