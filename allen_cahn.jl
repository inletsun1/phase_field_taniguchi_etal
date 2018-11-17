using PyPlot

#parameters
XSIZE = 50 + 2
YSIZE = 50 + 2
TSIZE = 100 
dx = 0.5e-6
dy = dx
gamma = 1.0
delta = 4*dx
M = 4.0e-14
rambda = 0.1
b = log((1-rambda)/rambda)

a = sqrt(3*delta*gamma/b)*10
W = 6*gamma*b/delta
Mp = b*M/(3*delta)
beta = 0.5

dt = dx^2/(5*Mp*a^2)
println("a^2 = ", a^2)
println("b = ", b)
println("W = ", W)
println("Mp = ", Mp)
println("dt = ", dt, "[s]")

#initialization
r0 = dx*XSIZE/2;
phi = zeros(Float64, YSIZE, XSIZE, TSIZE);
init_phi = [sqrt(i.^2 + j.^2) for i in range(0.0, stop=dx*XSIZE, length=XSIZE), j in range(0.0, stop=dy*YSIZE, length=YSIZE)] .- r0 
phi[:, :, 1] = (1.0 .- tanh.(sqrt(2.0 .* W)./(2.0 .*a).*init_phi))./2.0

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
Ay[1, 1:3] = [1.0 -2.0 1.0]
Ay[end, end-2:end] = [1.0 -2.0 1.0]
Ax[1:3, 1] = [1.0; -2.0; 1.0]
Ax[end-2:end, end] = [1.0; -2.0; 1.0]
#time development
for t in 2:TSIZE
    phi[:,:,t] = phi[:,:,t-1]  
    .+ dt .*Mp .*(4.0 .*W .*phi[:,:,t-1] .*(1 .-phi[:,:,t-1]) .*(phi[:,:,t-1] .-0.5 .+ beta) 
               .+ ((a ./dy) .^2) .* Ay * phi[:,:,t-1] 
               .+ ((a ./dx) .^2) .* phi[:,:,t-1] * Ax)
end

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

imshow(phi[:,:,1])
show()
imshow(phi[:,:,end])
show()
