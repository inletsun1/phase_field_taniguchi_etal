using PyPlot
using JLD


const XSIZE = 200
const YSIZE = 200
const TSIZE = 1000 
const dx = 0.1
const dy = dx
const Du = 0.05
const Dv = 0.20
const dt = 0.001

const S = 1.0
const tau = 0.83
const eps = 1.0
const eta = 1.0
const M = 0.5
const A0 = 5.0^2*pi
const R = 5.0

function main()
    #parameters
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

    function fu(U, V, phi)
        A = sum(phi.*dx.*dy)
        return phi.*(-alpha.*U.*V.^2 ./(Kk.+sum(phi.*V.^2 .*dx.*dy)./A) .+ beta.*U.*V./(Kp.+sum(U.*phi.*dx.*dy)./A) .+ S .-gamma.*U)
    end

    function fv(U, V, phi)
        A = sum(phi.*dx.*dy)
        return phi.*(alpha.*U.*V.^2 ./(Kk.+sum(phi.*V.^2 .*dx.*dy)./A) .- beta.*U.*V./(Kp.+sum(U.*phi.*dx.*dy)./A) .-mu.*V)
    end

    function fp(U, V, phi)
        return -eta./(tau.*eps.^2).*36.0.*phi.*(1.0.-phi).*(1.0.-2.0.*phi)
    end
    #initialization
    phi = zeros(Float64, YSIZE, XSIZE, TSIZE);
    U = zeros(Float64, YSIZE, XSIZE, TSIZE);
    V = zeros(Float64, YSIZE, XSIZE, TSIZE);
    phase_range = 0.01
    delta = dx*5
    k = log((1.0 - phase_range)/phase_range)/delta
    x0 = XSIZE/2.0 * dx
    y0 = YSIZE/2.0 * dy
    r0 = R + delta/2.0
    rmesh = [k*(sqrt((x*dx-x0)^2 + (y*dy-y0)^2)-r0) for y in 0:YSIZE-1, x in 0:XSIZE-1]
    phi[:,:,1] = (1.0 .- tanh.(rmesh))./2.0
    #V[:,:,1] = (1.0 .- tanh.(rmesh))./2.0
    #U[:,:,1] = 6.0.*(1.0 .- tanh.(rmesh))./2.0
    for x in 1:XSIZE
        for y in 1:YSIZE
            if ((x-1)*dx-x0)^2 + ((y-1)*dy-y0)^2<R-1.0
                U[y, x, 1] = 1.0/gamma
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
    for i in 2:YSIZE-1
        #By[i, i-1:i] = [-1.0 1.0] 
        By[i, i-1] = -1.0; By[i, i+1] = 1.0;
    end
    for i in 2:XSIZE-1
        #Bx[i-1:i, i] = [-1.0; 1.0] 
        Bx[i-1, i] = -1.0; Bx[i+1, i] = 1.0;
    end

    #Periodic
    Ay[1, 1:2] = [-2.0 1.0]; Ay[1, end] = 1.0;
    Ay[end, end-1:end] = [1.0 -2.0]; Ay[end, 1] = 1.0;
    Ax[1:2, 1] = [-2.0; 1.0]; Ax[end, 1] = 1.0;
    Ax[end-1:end, end] = [-2.0; 1.0]; Ax[1, end] = 1.0;
    #By[1, 1] = 1.0; By[1, end] = -1.0
    #Bx[1, 1] = 1.0; Bx[end, 1] = -1.0
    By[1,end] = -1.0; By[1,2] = 1.0; By[end, end-1] = -1.0; By[end, 1] = 1.0;
    By[end, 1] = -1.0; By[2,1] = 1.0; By[end-1, end] = -1.0; By[1, end] = 1.0;

    println("start calculation")
    #time development
    for t in 2:TSIZE
        A = sum(phi[:,:,t-1].*dx.*dy)

        abs_nabra_phi = sqrt.(((1.0 ./(2.0.*dx)) .* phi[:,:,t-1] *Bx).^2 .+ ((1.0 ./(2.0.*dy)) .* By *phi[:,:,t-1]).^2)

        del_phi = (eta.*((1.0 ./dy) .^2.0 .* Ay * phi[:,:,t-1] .+ (1.0 ./dx) .^2.0 .* phi[:,:,t-1] * Ax) .-M.*(A .- A0).*abs_nabra_phi .+ (a.*V[:,:,t-1] .- b.*U[:,:,t-1]).*abs_nabra_phi)./tau .+ fp(U[:,:,t-1], V[:,:,t-1], phi[:,:,t-1])

        phi[:,:,t] = phi[:,:,t-1] .+ dt.*del_phi

        reaction_U = fu(U[:,:,t-1], V[:,:,t-1], phi[:,:,t-1]).-Xu.*U[:,:,t-1].*abs_nabra_phi.^2 ./(sum(dx.*dy.*abs_nabra_phi.^2))

        reaction_V = fv(U[:,:,t-1], V[:,:,t-1], phi[:,:,t-1]).-Xv.*V[:,:,t-1].*abs_nabra_phi.^2 ./(sum(dx.*dy.*abs_nabra_phi.^2))

        U[:,:,t] = U[:,:,t-1].+ dt .*(Du.*(phi[:,:,t-1] .*((1.0 ./dy) .^2 .* Ay * U[:,:,t-1] .+ (1.0 ./dx) .^2 .* U[:,:,t-1] * Ax) .+ ((1.0 ./dx) .* phi[:,:,t-1] *Bx .+ (1.0 ./dy) .* By *phi[:,:,t-1]) .*((1.0 ./dx) .* U[:,:,t-1] *Bx .+ (1.0 ./dy) .* By *U[:,:,t-1])) .+ reaction_U .- del_phi)

        V[:,:,t] = V[:,:,t-1].+ dt .*(Dv.*(phi[:,:,t-1] .*((1.0 ./dy) .^2 .* Ay * V[:,:,t-1] .+ (1.0 ./dx) .^2 .* V[:,:,t-1] * Ax) .+ ((1.0 ./dx) .* phi[:,:,t-1] *Bx .+ (1.0 ./dy) .* By *phi[:,:,t-1]) .*((1.0 ./dx) .* V[:,:,t-1] *Bx .+ (1.0 ./dy) .* By *V[:,:,t-1])) .+ reaction_V .- del_phi)
    end
    println("end calculation")

    JLD.save("./simulation_results_rungekutta.jld", "phi", phi, "V", V, "U", U)
    println("finished to save the results")
end

main()
