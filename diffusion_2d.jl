using PyPlot

#parameters
XSIZE = 100 + 2
YSIZE = 100 + 2
TSIZE = 1000 
Lx = 100.0
Ly = 100.0
dx = Lx/(XSIZE-2)
dy = dx
D = 1
dt = 0.1
println("dt = ", dt, ", 1/(2D(1/dx^2 + 1/dy^2)) = ", 1/(2*D*(1/(dx^2) + 1/(dy^2))))


#initialization
r0 = dx*5;
phi = zeros(Float64, YSIZE, XSIZE, TSIZE);
for i in 1:YSIZE
    for j in 1:XSIZE
        if sqrt((i-YSIZE/2)^2 + (j-XSIZE/2)^2) <= r0/dx
            phi[i, j, 1] = 1.0
        end
    end
end

#boundary condition and diffusion discritization matrix
Ax = zeros(Float64, XSIZE, XSIZE)
Ay = zeros(Float64, YSIZE, YSIZE)
B = zeros(Float64, YSIZE, XSIZE)
for i in 2:YSIZE-1
    Ay[i, i-1:i+1] = [1.0 -2.0 1.0] 
end
for i in 2:XSIZE-1
    Ax[i-1:i+1, i] = [1.0; -2.0; 1.0] 
end
#Dirichlet

Ay[1, 1:3] = [1.0 -2.0 1.0]
Ay[end, end-2:end] = [1.0 -2.0 1.0]
Ax[1:3, 1] = [1.0; -2.0; 1.0]
Ax[end-2:end, end] = [1.0; -2.0; 1.0]

#=
#Periodic
Ay[1, 1:2] = [-2.0 1.0]; Ay[1, end] = 1.0;
Ay[end, end-1:end] = [1.0 -2.0]; Ay[end, 1] = 1.0;
Ax[1:2, 1] = [-2.0; 1.0]; Ax[end, 1] = 1.0;
Ax[end-1:end, end] = [-2.0; 1.0]; Ax[1, end] = 1.0;
=#
#time development

for t in 2:TSIZE
    phi[:,:,t] = phi[:,:,t-1].+ dt .* D.*((1 ./dy) .^2 .* Ay * phi[:,:,t-1] .+ (1 ./dx) .^2 .* phi[:,:,t-1] * Ax)
end

#imshow(phi[:, 51, :])
plot(phi[:,51,1])
plot(phi[:,51,Int(TSIZE/2)])
plot(phi[:,51,end])
show()
