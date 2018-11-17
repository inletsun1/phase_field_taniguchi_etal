using PyPlot

#parameters
XSIZE = 200 + 2
TSIZE = 500 
Lx = 100.0
dx = Lx/(XSIZE-2)
D = 1.0
dt = 0.1

println("dt = ", dt, ", dx^2/2D = ", dx^2/(2*D))

#initialization
phi = zeros(Float64, XSIZE, TSIZE);
#init_phi = [sqrt(i.^2 + j.^2) for i in range(0.0, stop=dx*XSIZE, length=XSIZE), j in range(0.0, stop=dy*YSIZE, length=YSIZE)] .- r0 
for i in (101-5):(101+5)
    phi[i, 1] = 1.0
end

#boundary condition and diffusion discritization matrix
Ax = zeros(Float64, XSIZE, XSIZE)
#Ay = zeros(Float64, YSIZE, YSIZE)
for i in 2:XSIZE-1
    Ax[i, i-1:i+1] = [1.0 -2.0 1.0] 
end
#=
for i in 2:XSIZE-1
    Ax[i-1:i+1, i] = [1.0; -2.0; 1.0] 
end
=#
#Dirichlet
#=
Ay[1, 1:3] = [1.0 -2.0 1.0]
Ay[end, end-2:end] = [1.0 -2.0 1.0]
Ax[1:3, 1] = [1.0; -2.0; 1.0]
Ax[end-2:end, end] = [1.0; -2.0; 1.0]
=#
#Periodic
Ax[1, 1:2] = [-2.0 1.0]; Ax[1, end] = 1.0;
Ax[end, end-1:end] = [1.0 -2.0]; Ax[end, 1] = 1.0;
#time development

for t in 2:TSIZE
    phi[:,t] = phi[:,t-1] .+ dt .* D.*(1 ./dx) .^2 .* Ax *phi[:,t-1]
end

#=
for t in 2:TSIZE
    for ny in 2:YSIZE-1
        for nx in 2:XSIZE-1
            phi[ny, nx, t] = phi[ny, nx, t-1] 
            + dt*D*(
                         (phi[ny-1, nx, t-1]-2*phi[ny, nx, t-1]+phi[ny+1, nx, t-1])/(dy^2) 
                         +(phi[ny, nx-1, t-1]-2*phi[ny, nx, t-1]+phi[ny, nx+1, t-1])/(dx^2))
        end
    end
    phi[1, :, t] = phi[2, :, t]
    phi[end, :, t] = phi[end-1, :, t]
    phi[:, 1, t] = phi[:, 2, t]
    phi[:, end, t] = phi[:, end-1, t]
end
=#
#=
println("Ax = ")
for i in 1:XSIZE
    println(Ax[i, :])
end
println("")
println("Ay = ")
for i in 1:YSIZE
    println(Ay[i, :])
end
=#
plot(phi[:,1])
plot(phi[:, Int(TSIZE/2)])
plot(phi[:, end])
show()
