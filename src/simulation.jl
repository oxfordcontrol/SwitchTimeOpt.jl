# Linear Systems STO
function simulate(m::linSTO)

  # Define Time
  t = collect(linspace(m.STOev.t0, m.STOev.tf, 10000))

  # Perform Simulation
  x, xpts, J = simulatelinsto(m.tau, m.STOev.x0, m.STOev.Q, m.STOev.E, m.STOev.A, t)

  return x, xpts, J, t

end

function simulate(m::linSTO, tau::Array{Float64,1})

  # Define Time
  t = collect(linspace(m.STOev.t0, m.STOev.tf, 10000))

  # Perform Simulation
  x, xpts, J = simulatelinsto(tau, m.STOev.x0, m.STOev.Q, m.STOev.E, m.STOev.A, t)

  return x, xpts, J, t

end


function simulate(m::linSTO, tau::Array{Float64,1}, t::Array{Float64,1})

  # Perform Simulation
  x, xpts, J = simulatelinsto(tau, m.STOev.x0, m.STOev.Q, m.STOev.E, m.STOev.A, t)

  return x, xpts, J, t

end






# Nonlinear Systems STO
function simulate(m::nlinSTO)

  # Define Time
  t = collect(linspace(m.STOev.t0, m.STOev.tf, 10000))

  # Get original Q, x0
  Q = m.STOev.Q[1:end-1, 1:end-1]
  E = m.STOev.E[1:end-1, 1:end-1]
  x0 = m.STOev.x0[1:end-1]

  # Perform Simulation
  x, xpts, J = simulatenlinsto(m.STOev.nonlin_dyn, m.tau, x0,  Q, E,  m.STOev.uvec, t)


  return x, xpts, J, t

end


function simulate(m::nlinSTO, tau::Array{Float64, 1})

  # Define Time
  t = collect(linspace(m.STOev.t0, m.STOev.tf, 10000))

  # Get original, Q, x0
  Q = m.STOev.Q[1:end-1, 1:end-1]
  E = m.STOev.E[1:end-1, 1:end-1]
  x0 = m.STOev.x0[1:end-1]

  # Perform Simulation
  x, xpts, J = simulatenlinsto(m.STOev.nonlin_dyn, tau, x0,  Q, E,  m.STOev.uvec, t)


  return x, xpts, J, t

end


function simulate(m::nlinSTO, tau::Array{Float64, 1}, t::Array{Float64, 1})

  # Get original, Q, x0
  Q = m.STOev.Q[1:end-1, 1:end-1]
  E = m.STOev.E[1:end-1, 1:end-1]
  x0 = m.STOev.x0[1:end-1]

  # Perform Simulation
  x, xpts, J = simulatenlinsto(m.STOev.nonlin_dyn, tau, x0,  Q, E,  m.STOev.uvec, t)

  return x, xpts, J, t

end


# Linearized Nonlinear System STO
function simulatelinearized(m::nlinSTO)

  # Define Time
  t = collect(linspace(m.STOev.t0, m.STOev.tf, 10000))

  # Get original Q, x0
  Q = m.STOev.Q[1:end-1, 1:end-1]
  E = m.STOev.E[1:end-1, 1:end-1]
  x0 = m.STOev.x0[1:end-1]

  # Perform Simulation
  x, xpts, J = simulatelinearizedsto(m.STOev.nonlin_dyn, m.STOev.nonlin_dyn_deriv, m.tau, m.STOev.tgrid, m.STOev.tfdelta, m.STOev.uvec,  x0, Q, E,  t)

  return x, xpts, J, t

end


# Linearized Nonlinear System STO
function simulatelinearized(m::nlinSTO, tau::Array{Float64,1})

  # Define Time
  t = collect(linspace(m.STOev.t0, m.STOev.tf, 10000))

  # Get original Q, x0
  Q = m.STOev.Q[1:end-1, 1:end-1]
  E = m.STOev.E[1:end-1, 1:end-1]
  x0 = m.STOev.x0[1:end-1]

  # Perform Simulation
  x, xpts, J = simulatelinearizedsto(m.STOev.nonlin_dyn, m.STOev.nonlin_dyn_deriv, tau, m.STOev.tgrid, m.STOev.uvec,  x0, Q, E,  t)

  return x, xpts, J, t

end




# Linearized Nonlinear System STO
function simulatelinearized(m::nlinSTO, tau::Array{Float64,1}, t::Array{Float64, 1})

  # Get original Q, x0
  Q = m.STOev.Q[1:end-1, 1:end-1]
  E = m.STOev.E[1:end-1, 1:end-1]
  x0 = m.STOev.x0[1:end-1]

  # Perform Simulation
  x, xpts, J = simulatelinearizedsto(m.STOev.nonlin_dyn, m.STOev.nonlin_dyn_deriv, tau, m.STOev.tgrid, m.STOev.uvec, x0, Q, E, t)

  return x, xpts, J, t

end


### Lower Level Functions

# Linear Systems
function simulatelinsto(tau::Array{Float64,1}, x0::Array{Float64, 1}, Q::Array{Float64,2}, E::Array{Float64,2}, A::Array{Float64,3}, t::Array{Float64,1})
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
  Jtoint = diag(x'*Q*x)
  J = trapz(t, Jtoint) + (x[:, end]'*E*x[:, end])[1]

  return x, xpts, J


end


# Linearized Nonlinear System
  function simulatelinearizedsto(nonlin_dyn::Function, nonlin_dyn_deriv::Function, tau::Array{Float64,1}, tgrid::Array{Float64,1}, uvec::Array{Float64, 2}, x0::Array{Float64, 1}, Q::Array{Float64,2},  E::Array{Float64,2}, t::Array{Float64,1})

    # Get dimensions
    nx = length(x0)  # Number of States
    N = length(tau)  # Number of switches
    ngrid = length(tgrid)  # Number of elements in time grid

    # Create merged and sorted time vector with grid and switching times
    tvec, tauIdx = mergeSortFindIndex(tgrid, tau)

    # Create matrix of Linearized Dynamics
    A = Array(Float64, nx+1, nx+1, N+ngrid - 1)

    # Compute Initial States
    xpts = Array(Float64, nx+1, N+ngrid)  # Augmented State for Linearization
    xpts[:,1] = [x0; 1]


      uIdx = 1  # Initialize index for current u

      # Compute Matrix Exponentials
      for i =1:N+ngrid-1  # Iterate over all grid (incl sw times)

        # Verify which U input applies
        if uIdx <= N
          if i>= tauIdx[uIdx + 1]
            uIdx += 1
          end
        end

        # Linearize Dynamics
        A[:,:,i] = linearizeDyn(nonlin_dyn, nonlin_dyn_deriv, xpts[1:end-1,i], uvec[:,uIdx])

        # Compute Next Point in Simulation
        xpts[:,i+1] = expm(A[:, :, i]*(tvec[i+1] - tvec[i]))*xpts[:, i]

      end


      # Compute State Trajectory
      x = Array(Float64, nx+1, length(t))
      x[:, 1] = [x0; 1]
      tauind = 1  # Index to keep track of the current mode


      for i = 2:length(t)

        # Check if we are still in the current switching interval. Otherwise Change
        if tauind < N+ngrid-1
          if t[i] > tvec[tauind + 1]
            tauind += 1
          end
        end

        # Compute State
        x[:, i] = expm(A[:, :, tauind]*(t[i] - tvec[tauind]))*xpts[:, tauind]

      end


      # Remove Extra State for Cost Function Evaluation
      x = x[1:end-1, :]

      # Numerically Integrate Cost Function
      Jtoint = diag(x'*Q*x)
      J = trapz(t, Jtoint) + (x[:, end]'*E*x[:, end])[1]

      return x, xpts, J

  end



# Noninear Systems
function simulatenlinsto(nonlin_dyn::Function, tau::Array{Float64,1}, x0::Array{Float64, 1}, Q::Array{Float64,2}, E::Array{Float64,2}, uvec::Array{Float64,2}, t::Array{Float64,1})
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

        _, xmap = ode45(nldyn, xprevSwitch, [t[tempInd1]; (t[tempInd1]+t[tempInd2])/2; t[tempInd2]], points=:specified)  # Integrate
        xmap = hcat(xmap...)
        xtemp = [xmap[:,1] xmap[:,3]]          # Take only two points
      else

        _, xmap = ode45(nldyn, xprevSwitch, t[tempInd1:tempInd2], points=:specified)
        xtemp = hcat(xmap...)  # https://github.com/JuliaLang/ODE.jl/issues/80

      end

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
  Jtoint = diag(x'*Q*x)
  J = trapz(t, Jtoint) + (x[:,end]'*E*x[:,end])[1]

  return x, xpts, J


end


#### Auxiliary functions
function simulateinput(m::nlinSTO)

  # Define Time
  t = collect(linspace(m.STOev.t0, m.STOev.tf, 10000))


  # u = computenlswinput(m.taucomplete, m.STOev.uvec, t)
  u = computenlswinput(m.tau, m.STOev.uvec, t)  # Fix all TAUCOMPLETE

  return u, t

end

function simulateinput(m::nlinSTO, t::Array{Float64, 1})

  # u = computenlswinput(m.taucomplete, m.STOev.uvec, t)
  u = computenlswinput(m.tau, m.STOev.uvec, t)  # Fix all TAUCOMPLETE


  return u, t

end



# Compute Actual Inputs from Artificial Ones
function computenlswinput(tauopt, uvec, t)
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
    return r/2
end
