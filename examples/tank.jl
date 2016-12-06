
using SwitchTimeOpt

# Plotting
using PyPlot, PyCall
close("all")
# Import Seaborn for nice plots
@pyimport seaborn as sns
sns.set_palette("hls", 8)
sns.set_context("paper", font_scale=1.5)
sns.set_style("whitegrid")
# Use Latex Labels in Plots
plt[:rc]("text", usetex=true)
plt[:rc]("font", family="serif")   # Use Serif math


maxiter = 15
using Ipopt
solver = IpoptSolver(
  print_level=0,
  tol = 1e-10,  # Increase required precision to avoid convergence before maxiter iterations
  max_iter = maxiter,
  linear_solver="ma57")


### Define system parameters
#  Time vector
t0 = 0.0
tf = 10.0

# Integer inputs
uvec = repmat([1.0; 2.0], 8, 1)'

# Number of switches
N = size(uvec,2) -1

# Cost function matrix
K = 1.0
C = sqrt(K)*[0.0 1.0 -1.0]
Q = C'*C

# Initial state
x0 = [2.0; 2.0; 3.0]


### Define system dynamics
function nldyn(x::Array{Float64,1}, u::Array{Float64,1})
  n = length(x)
  f = zeros(n)

  f[1] = -sqrt(x[1]) + u[1]
  f[2] = sqrt(x[1]) - sqrt(x[2])
  f[3] = -0.05
  return f

end

function nldyn_deriv(x::Array{Float64,1}, u::Array{Float64,1})
  df = [-1/(2*sqrt(x[1]))        0                      0;
        1/(2*sqrt(x[1]))        -1/(2*sqrt(x[2]))       0;
        0                        0                      0]
end




# Generate and solve problems with different grid points
ngrid = [10; 30; 50; 100]


# Preallocate vectors for results
objode45 = Array(Float64, length(ngrid))
objlin = Array(Float64, length(ngrid))
objiterates = Array(Float64, maxiter+1, length(ngrid))
cputime = Array(Float64, length(ngrid))
nobjeval = Array(Int, length(ngrid))
ngradeval = Array(Int, length(ngrid))
nhesseval = Array(Int, length(ngrid))
tauopt = Array(Float64, N, length(ngrid))
xsim = Array(Float64, 3, 10^4, length(ngrid))
xlinsim = Array(Float64, 3, 10^4, length(ngrid))
usim = Array(Float64, 1, 10^4, length(ngrid))

# Initialize the problem first
m = stoproblem(x0, nldyn, nldyn_deriv, uvec)


for i = 1:length(ngrid) # Iterate over all grid points numbers
  m = stoproblem(
    x0,                 # Initial state
    nldyn,              # Nonlinear dynamics
    nldyn_deriv,        # Nonlinear dynamics derivative
    uvec,               # Vector of integer inputs
    ngrid = ngrid[i],   # Number of grid points
    t0 = t0,            # InitialtTime
    tf = tf,            # Final time
    Q = Q,              # Cost matrix
    solver = solver)


  # Solve optimization
  solve!(m)

  # Obtain values
  tauopt[:, i] = gettau(m)
  objlin[i] = getobjval(m)
  cputime[i] = getsoltime(m)
  nobjeval[i] = getnobjeval(m)
  ngradeval[i] = getngradeval(m)
  nhesseval[i] = getnhesseval(m)


  # Simulate System
  # Nonlinear Simulation
  xsim[:, :, i], xpts, objode45[i], t = simulate(m, tauopt[:, i])
  usim[:, :, i], _ = simulateinput(m, t)

  # Linearized simulation
  xlinsim[:, :, i], _, Jlinsim, _ = simulatelinearized(m, tauopt[:, i], t)

  # Save objective iterates
  objiterates[:, i] = m.STOev.obj[2:end]

end



###  Print results
@printf("RESULTS\n")
@printf("-------\n\n")

@printf("+==================================================================================+\n")
@printf("|  ngrid   |   objode45  |   objlin    | deltaobj [%%] |  nobjeval  |  cputime [s]  |\n")
@printf("+==================================================================================+\n")


for i = 1:length(ngrid)
  @printf("|  %7.i | %9.4f   | %9.4f   |   %6.3f     | %9.i  | %12.4f  | \n", ngrid[i], objode45[i], objlin[i], 100*norm(objode45[i]- objlin[i])/objode45[i], nobjeval[i], cputime[i])
  @printf("+----------------------------------------------------------------------------------+\n")
end



@printf("\nOptimal Switching Times for ngrid = 10\n")
@printf("--------------------------------------\n")
@printf("tauopt = "); show(round(tauopt[:, 1],2)); @printf("\n")



# Generate plots for ngrid = 10
t = linspace(t0, tf, 10000)

figure()
subplot(3,1,1)
plot(t, xlinsim[1,:, 1], sns.xkcd_rgb["grass green"], linestyle = "dashdot")
plot(t, xsim[1,:, 1], sns.xkcd_rgb["denim blue"])
ylim(1.5, 4)
yticks([2; 3])
ylabel(L"x_1")


subplot(3,1,2)
plot(t, xsim[3,:, 1], sns.xkcd_rgb["black"], linestyle = "dotted")
plot(t, xlinsim[2,:, 1], sns.xkcd_rgb["grass green"], linestyle = "dashdot")
plot(t, xsim[2,:, 1], sns.xkcd_rgb["denim blue"])
ylabel(L"x_2")
ylim(1.5, 4.0)
yticks([2; 3])

subplot(3,1,3)
ax1 = plot(t, usim[1,:, 1], sns.xkcd_rgb["denim blue"])
ylim(0.8, 2.2)
yticks([1; 2])
ax = gca()
ax[:set_yticklabels]([L"u_{\mathrm{min}}", L"u_{\mathrm{max}}"])
ylabel(L"u")
xlabel(L"$\mathrm{Time}\; [s]$")

## Save figure
# tight_layout()
# savefig("tank_problem.pdf")


# Do not return anything
nothing
