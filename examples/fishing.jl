using SwitchTimeOpt

# Plotting settings
using PyPlot, PyCall
close("all")
# Import seaborn for nice plots
@pyimport seaborn as sns
sns.set_palette("hls", 8)
sns.set_context("paper", font_scale=1.5)
sns.set_style("whitegrid")
# Use Latex Labels in Plots
plt[:rc]("text", usetex=true)
plt[:rc]("font", family="serif")

maxiter = 25
using Ipopt
solver = IpoptSolver(
  print_level=0,  # Suppress output
  # alpha_red_factor = 0.05,  # Reduction factor for line search step
  max_iter = maxiter,
  linear_solver="ma57")
# using KNITRO
# solver = KnitroSolver(KTR_PARAM_MAXIT = 25,
# # KTR_PARAM_FEASTOL=1e-02,
# # KTR_PARAM_INFEASTOL=1e-2,
# # KTR_PARAM_OPTTOL=1e-03,
# # KTR_PARAM_BAR_FEASMODETOL = 1e-03,
# # KTR_PARAM_DERIVCHECK_TOL=1e-03
# # KTR_PARAM_FEASTOLABS=1e-02,
# # KTR_PARAM_MSSAVETOL=1e-4,
# KTR_PARAM_ALG = 4
# )

### Define system parameters
# Time vector
t0 = 0.0
tf = 12.0

# Size of the system
nx = 4

# Integer input
uvec = [repmat([0.0; 1.0], 4, 1); 0.0]'

# Number of switching times
N = size(uvec,2) - 1

# Cost funcction matrix
C = [1.0 0.0 -1.0 0.0;
     0.0 1.0 0.0  -1.0]
Q = C'*C

# Define initial state
x0 = [0.5; 0.7; 1; 1]


### Define system dynamics
function nldyn(x::Array{Float64,1}, u::Array{Float64,1})
  n = length(x)
  f = zeros(n)
  f[1] = x[1] - x[1]*x[2] - 0.4*x[1]*u[1]
  f[2] = -x[2] + x[1]*x[2] - 0.2*x[2]*u[1]
  return f

end

function nldyn_deriv(x::Array{Float64,1}, u::Array{Float64,1})
  df = [1.0-x[2]-0.4*u[1]       -x[1]                   0    0;
        x[2]                     -1+x[1]-0.2*u[1]       0    0;
        0                        0                      0    0;
        0                        0                      0    0]
end


### Generate and solve problems with different grid points
ngrid = [25; 50; 100; 150; 200; 250; 300]

# Preallocate vectors for results
objode45 = Array(Float64, length(ngrid))
objlin = Array(Float64, length(ngrid))
objiterates = Array(Float64, maxiter+1, length(ngrid))
cputime = Array(Float64, length(ngrid))
nobjeval = Array(Int, length(ngrid))
ngradeval = Array(Int, length(ngrid))
nhesseval = Array(Int, length(ngrid))
tauopt = Array(Float64, N, length(ngrid))
xsim = Array(Float64, 4, 10^4, length(ngrid))
xlinsim = Array(Float64, 4, 10^4, length(ngrid))
usim = Array(Float64, 1, 10^4, length(ngrid))

# Initialize the problem first
m = stoproblem(x0, nldyn, nldyn_deriv, uvec)

for i = 1:length(ngrid)  # Iterate over all grid points numbers

  m = stoproblem(
    x0,                 # Initial state
    nldyn,              # Nonlinear dynamics
    nldyn_deriv,        # Nonlinear dynamics derivative
    uvec,               # Vector of integer inputs
    ngrid = ngrid[i],   # Number of grid points
    t0 = t0,            # Initial time
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

  # Simulate system
  # Nonlinear simulation
  xsim[:, :, i], xpts, objode45[i], t = simulate(m, tauopt[:, i])
  usim[:, :, i], _ = simulateinput(m, t)

  # Linearized simulation
  xlinsim[:, :, i], _, Jlinsim, _ = simulatelinearized(m, tauopt[:, i], t)

  # Save objective function iterates
  objiterates[:, i] = m.STOev.obj[2:end]

end


###  Print results
@printf("RESULTS\n")
@printf("-------\n\n")

@printf("+=================================================================================+\n")
@printf("|  ngrid   |   objode45  |   objlin    |   deltaobj  |  nobjeval  |  cputime [s]  |\n")
@printf("+=================================================================================+\n")


for i = 1:length(ngrid)
  @printf("|  %7.i | %9.4f   | %9.4f   | %6.2fe-03  | %9.i  | %12.4f  | \n", ngrid[i], objode45[i], objlin[i], 10^3*norm(objode45[i]- objlin[i]), nobjeval[i], cputime[i])
  @printf("+---------------------------------------------------------------------------------+\n")
end



@printf("\nOptimal Switching Times for ngrid = 25\n")
@printf("--------------------------------------\n")
@printf("tauopt = "); show(round(tauopt[:, 1],3)); @printf("\n")


# # Generate table content for latex file
# Mtowrite = [ngrid objode45 objlin nobjeval cputime]
# f = open("Mfishproblem.csv","w")
# for i = 1:length(ngrid)
#   @printf(f, "%i & %.4f & %.4f & %.2fe-03 & %i & %.2f\\\\\n", ngrid[i], objode45[i], objlin[i],  10^3*norm(objode45[i]- objlin[i]), nobjeval[i], cputime[i])
# end
# close(f)


# Generate plots for ngrid = 25
t = linspace(t0, tf, 10000)

figure()
subplot(3,1,1)
plot(t, xlinsim[1,:, 1]', sns.xkcd_rgb["grass green"], linestyle = "dashdot")
plot(t, xsim[1,:, 1]', sns.xkcd_rgb["denim blue"])
ylim(0, 1.75)
yticks([0; 1; ])
ylabel(L"x_1")

subplot(3,1,2)
plot(t, xlinsim[2,:, 1]', sns.xkcd_rgb["grass green"], linestyle = "dashdot")
plot(t, xsim[2, :, 1]', sns.xkcd_rgb["denim blue"])
ylabel(L"x_2")
ylim(0, 1.75)
yticks([0; 1; ])

subplot(3,1,3)
plot(t, usim[1,:, 1]', sns.xkcd_rgb["denim blue"])
ylim(-0.2, 1.2)
yticks([0; 1])
ylabel(L"u")
xlabel(L"$\mathrm{Time}\; [s]$")

# Save figure
tight_layout()
savefig("fishing_problem.pdf")



# Plot cost function iterates for all the grid numbers
figure()
for i = 1:length(ngrid)
  stemp = latexstring(@sprintf("n_{\\text{grid}} = %i", ngrid[i]))
  semilogy(objiterates[:, i], label=stemp)
end
legend()
yticks([1; 2; 3; 4; 5; 6])
ylabel(L"$J$")
xlabel(L"\# iterations")
tight_layout()
savefig("Jfishing_problem.pdf")


# Do not return anything
nothing
