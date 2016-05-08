# Trapezoidal integration rule
function trapz{Tx<:Number, Ty<:Number}(x::Vector{Tx}, y::Vector{Ty})
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end
    r = zero(zero(Tx) + zero(Ty))
    if n == 1; return r; end
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    #= correction -h^2/12 * (f'(b) - f'(a))
    ha = x[2] - x[1]
    he = x[end] - x[end-1]
    ra = (y[2] - y[1]) / ha
    re = (y[end] - y[end-1]) / he
    r/2 - ha*he/12 * (re - ra)
    =#
    return r/2
end

# Linear Systems
function simulate(m::linSTO)

  # Define Time
  t = collect(linspace(m.STOev.t0, m.STOev.tf, 10000))

  # Perform Simulation
  x, xpts, J = simulateLinSTO(m.tau, m.STOev.x0, m.STOev.Q, m.STOev.A, t)

  return x, xpts, J, t

end

# Nonlinear Systems
function simulate(m::nlinSTO)

  # Define Time
  t = collect(linspace(m.STOev.t0, m.STOev.tf, 10000))

  # Get original uvec, Q, x0
  uvec = m.STOev.uvec[:,m.nartsw:m.nartsw:end]
  Q = m.STOev.Q[1:end-1, 1:end-1]
  x0 = m.STOev.x0[1:end-1]

  # Define inputs in time
  u = computeNlSwInput(m.tau, uvec, t);

  # Perform Simulation
  x, xpts, J = simulateNlinSTO(m.STOev.nonlin_dyn, m.tau, x0,  Q,  uvec, t)


  return x, xpts, u, J, t

end

# Linearized Nonlinear System
function simulatelinearized(m::nlinSTO)

  # Define Time
  t = collect(linspace(m.STOev.t0, m.STOev.tf, 10000))

  # Define inputs in time
  u = computeNlSwInput(m.taucomplete, m.STOev.uvec, t);

  # Get original Q, x0
  Q = m.STOev.Q[1:end-1, 1:end-1]
  x0 = m.STOev.x0[1:end-1]

  # Perform Simulation
  x, xpts, J = simulateLinearizedSTO(m.STOev.nonlin_dyn, m.STOev.nonlin_dyn_deriv, m.taucomplete, m.STOev.uvec,  x0, Q,  t)

  return x, xpts, u, J, t

end


function simulateLinSTO(tau::Array{Float64,1}, x0::Array{Float64, 1}, Q::Array{Float64,2}, A::Array{Float64,3}, t::Array{Float64,1})
  # Get dimensions
  nx = size(A, 1)  # Number of States
  N = length(tau)  # Number of switches

  tau = [t[1]; tau; t[end]]  # Extend tau vector to simplify numbering

  # Compute Initial States
  xpts = Array(Float64, nx, N+1)
  xpts[:,1] = x0
  for i = 2:N+1
    xpts[:,i] = expm(A[:, :, i-1]*(tau[i] - tau[i-1]))*xpts[:, i-1]
  end

  # Compute State Trajectory
  x = Array(Float64, nx, length(t))
  x[:, 1] = x0
  tauind = 1  # Index to keep track of the current mode

  for i = 2:length(t)
    # Check if we are still in the current switching interval. Otherwise Change
    if tauind < N+1
      if t[i] > tau[tauind + 1]
        tauind += 1
      end
    end

    # Compute State
    x[:, i] = expm(A[:, :, tauind]*(t[i] - tau[tauind]))*xpts[:, tauind]

  end

  # Numerically Integrate Cost Function
  Jtoint = 1/2*diag(x'*Q*x)
  J = trapz(t, Jtoint)

  return x, xpts, J


end


# Linearized Nonlinear System
  function simulateLinearizedSTO(nonlin_dyn::Function, nonlin_dyn_deriv::Function, tau::Array{Float64,1}, uvec::Array{Float64, 2}, x0::Array{Float64, 1}, Q::Array{Float64,2},  t::Array{Float64,1})

    # # Check if nonlinear dynamics and derivatives are defined
    # try
    #   method_exists(nonlin_dyn, (Array{Float64}, Array{Float64}, Array{Float64}))
    # catch err
    #   println("No function nonlin_dyn defined. You need to define the dynamics of your system to perform switching time optimization")
    #   throw(err)
    # end
    #
    # try
    #   method_exists(nonlin_dyn_deriv, (Array{Float64}, Array{Float64}, Array{Float64}))
    # catch err
    #   println("No function nonlin_dyn_deriv defined. You need to define the jacobian of the dynamics of your system with respect to the state variable to perform switching time optimization")
    #   throw(err)
    # end

    # Get dimensions
    nx = length(x0)  # Number of States
    N = length(tau)  # Number of switches

    tau = [t[1]; tau; t[end]]  # Extend tau vector to simplify numbering

    # Create matrix of Linearized Dynamics
    A = Array(Float64, nx+1, nx+1, N+1)

    # Compute Initial States
    xpts = Array(Float64, nx+1, N+1)  # Augmented State for Linearization
    xpts[:,1] = [x0; 1]

    for i = 2:N+1

      # Generate Linearized Dynamics
      A[:,:,i-1] = linearizeDyn(nonlin_dyn, nonlin_dyn_deriv, xpts[1:end-1,i-1], uvec[:,i-1])

      # Compute Next Point in Simulation
      xpts[:,i] = expm(A[:, :, i-1]*(tau[i] - tau[i-1]))*xpts[:, i-1]

    end

    # Generate Linearized Dynamics for last input
    A[:,:,N+1] = linearizeDyn(nonlin_dyn, nonlin_dyn_deriv, xpts[1:end-1,N+1], uvec[:,N+1])

    # Compute State Trajectory
    x = Array(Float64, nx+1, length(t))
    x[:, 1] = [x0; 1]
    tauind = 1  # Index to keep track of the current mode

    for i = 2:length(t)
      # Check if we are still in the current switching interval. Otherwise Change
      if tauind < N+1
        if t[i] > tau[tauind + 1]
          tauind += 1
        end
      end

      # Compute State
      x[:, i] = expm(A[:, :, tauind]*(t[i] - tau[tauind]))*xpts[:, tauind]

    end

    # Remove Extra State for Cost Function Evaluation
    x = x[1:end-1, :]

    # Numerically Integrate Cost Function
    Jtoint = 1/2*diag(x'*Q*x)
    J = trapz(t, Jtoint)

    return x, xpts, J


  end



# Noninear Systems
function simulateNlinSTO(nonlin_dyn::Function, tau::Array{Float64,1}, x0::Array{Float64, 1}, Q::Array{Float64,2}, uvec::Array{Float64,2}, t::Array{Float64,1})
  # Get dimensions
  nx = length(x0)  # Number of States
  N = length(tau)  # Number of switches

  tau = [t[1]; tau; t[end]]  # Extend tau vector to simplify numbering

  # Define empty state trajectory
  x = zeros(nx, length(t))
  x[:,1] = x0

  # Define indeces to determine current switching mode
  tempInd1 = 1
  tempInd2 = 1
  xprevSwitch = x0

  # Create Vector of Points at the switching instants
  xpts = zeros(nx, N+1)
  xpts[:, 1] = x0

  for i = 1:N+1 # Integrate over all the intervals

    # redefine Dynamic function
    nldyn(t, x) = nonlin_dyn(x, uvec[:, i])

    while t[tempInd2] < tau[i+1]
      tempInd2 = tempInd2 + 1  # Increase time index
    end

    if tempInd2>tempInd1  # There has been a progress and the switching times are not collapsed. So we integrate.
      if tempInd2 == tempInd1 + 1  # Only one step progress. Handle as special case for the integrator
        # println("$(tempInd1)")
        # println("$(tempInd2)")
        # show([t[tempInd1] (t[tempInd1]+t[tempInd2])/2 t[tempInd2]])
        _, xmap = ode45(nldyn, xprevSwitch, [t[tempInd1]; (t[tempInd1]+t[tempInd2])/2; t[tempInd2]], points=:specified)  # Integrate
        xmap = hcat(xmap...)
        xtemp = [xmap[:,1] xmap[:,3]]          # Take only two points
      else
        # println("$(tempInd1)")
        # println("$(tempInd2)")
        # show(t[tempInd1:tempInd2])
        _, xmap = ode45(nldyn, xprevSwitch, t[tempInd1:tempInd2], points=:specified)
        xtemp = hcat(xmap...)  # https://github.com/JuliaLang/ODE.jl/issues/80
        # xtemp = xmap[1]'
        # for j = 2:nx
        #   xtemp = [xtemp; xmap[j]']
        # end
      end

      # http://stackoverflow.com/questions/24275980/slice-array-of-arrays-in-julia
      x[:, tempInd1:tempInd2] = xtemp

      # Update indeces for next iteration
      xprevSwitch = x[:, tempInd2]
      tempInd1 = tempInd2

      # Update vector of switching states
      if i < N+1
        xpts[:, i+1] = x[:, tempInd2]
      end
    end


  end

  # Numerically Integrate Cost Function
  Jtoint = 1/2*diag(x'*Q*x)
  J = trapz(t, Jtoint)

  return x, xpts, J


end



# Compute Actual Inputs from Artificial Ones
function computeNlSwInput(tauopt, uvec, t)
  N = length(tauopt)
  nu = size(uvec, 1)
  tau = [t[1]; tauopt; t[end]]  # Add initial and final switching to simplify numbering

  # Define empty input vector
  usim = zeros(nu, length(t))

  # Define indeces to span time vector
  tempInd1 = 1
  tempInd2 = 1

  for i = 1:N+1  # Iterate over all intervals
    while t[tempInd2] < tau[i+1]
      tempInd2 += 1
    end

    if tempInd2 > tempInd1
      usim[:, tempInd1:tempInd2] = repmat(uvec[:, i], 1, tempInd2 - tempInd1 + 1)
      tempInd1 = tempInd2
    end


  end
  return usim


end
