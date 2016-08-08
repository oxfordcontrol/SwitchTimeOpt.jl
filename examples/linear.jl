using SwitchTimeOpt

# Plotting
using PyPlot, PyCall
close("all")
# Import Seaborn for nice plots
@pyimport seaborn as sns
sns.set_context("paper", font_scale=1.5)  # Set default style seaborn
sns.set_style("whitegrid")
# Use Latex Labels in Plots
plt[:rc]("text", usetex=true)
plt[:rc]("font", family="serif")   # Use Serif math

using Ipopt
solver = IpoptSolver(
          print_level = 0,
          tol = 1e-08,
          linear_solver="ma57")


### Define system parameters
# Size of the state space
nx = 2;

# Number of switches
N = 5;

# Cost function matrix
Q = 1*eye(2)
Qf = 0*eye(2)

# Initial state
x0 = [1.0; 1.0]

# Define initial and final time
t0 = 0.0
tf = 1.0


### Define system dynamics
A = zeros(nx, nx, N+1)

A[:, :, 1] = [-1 0;
              1  2]
A[:, :, 2] = [1 1;
              1 -2]
for i = 3:N+1
  A[:, :, i] = A[:, :, mod(i+1, 2)+1]
end

### Define and solve switching time optimization problem
m = stoproblem(x0,
              A,
              Q = Q,
              solver=solver)

# Solve problem
solve!(m)

# Obtain values
tauopt = gettau(m)
Jopt = getobjval(m)
soltime = getsoltime(m)

# Simulate linear system
xsim, xpts, Jsim, t = simulate(m)


### Print results
@printf("OBTAINED RESULTS\n")
@printf("----------------\n")
@printf("Objective Function at the Optimum:              J = %.3f\n", Jopt)
@printf("Simulated Objective Function at the Optimum: Jsim = %.3f\n", Jsim)
@printf("Elapsed Time:                                time = %.2f [ms]\n", mean(timeVec)*10^3)
@printf("Optimum:                                      tau = "); show(round(tauopt,3)); @printf("\n")

### Plot results
figure()
subplot(2,1,1)
for i = 1:m.STOev.ngrid
  axvline(x=m.STOev.tgrid[i], color="lightgrey", linestyle="--")
end
plot(t, xsim[1,:]')
ylabel(L"x_1")
yticks([1; 1.5; 2; 2.5])
ylim(0.5, 3)

subplot(2,1,2)
for i = 1:m.STOev.ngrid
  axvline(x=m.STOev.tgrid[i], color="lightgrey", linestyle="--")
end
plot(t, xsim[2,:]')
yticks([1; 1.5; 2; 2.5])
ylim(0.5, 3)
ylabel(L"x_2")
xlabel(L"$\mathrm{Time}\; [s]$")

tight_layout()
savefig("linearExample.pdf")


# Plot cost function iterates
figure()
plot(m.STOev.obj)
ylabel(L"$J$")
xlabel(L"\# iterations")




nothing  # Don't print anything
