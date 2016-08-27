# Define initialization Function
function MathProgBase.initialize(d::STOev, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in [:Grad, :Jac, :Hess])
            error("Unsupported feature $feat")
        end
    end
end


# List Available Features
MathProgBase.features_available(d::STOev) = [:Grad, :Jac, :Hess]




# Functions to precompute matrices for linear STO

"Precompute matrices for cost function"
function precompMatJ!(d::STOev, x)

  # Propagate dynamics
  propagateDynamics!(d, x)

  # DEBUG: Timing
  # timeprop = @elapsed propagateDynamics!(d, x)
  # @printf("Time elapsed propag = %.4f mus\n", timeprop*10^6)

  # Compute Matrices S
  computeSlw!(d)

  # DEBUG: Timing
  # timeS = @elapsed computeSlw!(d)
  # @printf("Time elapsed S = %.4f mus\n", timeS*10^6)

end

"Precompute matrices for gradient"
function precompMatGradJ!(d::STOev, x)
  if d.deltafun_prev!= x  # Different point than the previous cost function evaluation
    precompMatJ!(d, x)  # precompute matrices for cost function
    d.deltafun_prev[:] = x

    computeClw!(d)  # precompute matrices for gradient

    # DEBUG: Timing
    # timeC = @elapsed computeClw!(d)  # precompute matrices for gradient
    # @printf("Time elapsed C = %.4f mus\n", timeC*10^6)

  else
    computeClw!(d)  # precompute matrices for gradient

    # DEBUG: Timing
    # timeC = @elapsed computeClw!(d)  # precompute matrices for gradient
    # @printf("Time elapsed C = %.4f mus\n", timeC*10^6)
  end
end


"Precompute matrices for hessian"
function precompMatHessJ!(d::STOev, x)
  if d.deltagrad_prev!= x  # Different point than the previous gradient evaluation
    precompMatGradJ!(d, x)  # precompute matrices for gradient function
    d.deltagrad_prev[:] = x

    computePhilw!(d)  # precompute matrices for hessian

    # DEBUG: Timing
    # timePhi = @elapsed computePhilw!(d)  # precompute matrices for hessian
    # @printf("Time elapsed Phi = %.4f mus\n", timePhi*10^6)

  else
    computePhilw!(d)  # precompute matrices for hessian

    # DEBUG: Timing
    # timePhi = @elapsed computePhilw!(d)  # precompute matrices for hessian
    # @printf("Time elapsed Phi = %.4f mus\n", timePhi*10^6)

  end
end







# Evaluate Cost Function
function MathProgBase.eval_f(d::STOev, x)
  # Increase counter of objective function evaluations
  d.nobjeval += 1

  # Check if the matrices have already been precomputed
  if d.deltafun_prev!= x
    precompMatJ!(d, x)
    d.deltafun_prev[:] = x
  end

  J = (d.x0'*d.S[:,:,1]*d.x0)[1]

  return J
end



# Evaluate Gradient
function MathProgBase.eval_grad_f(d::STOev, grad_f, x)

  # Increase counter of gradient evaluations
  d.ngradeval += 1

  # Check if the matrices have already been precomputed
  if d.deltagrad_prev!= x
    precompMatGradJ!(d, x)
    d.deltagrad_prev[:] = x
  end


  # Update objective and cost function values to check iteration status (DEBUG)
  d.obj = [d.obj; MathProgBase.eval_f(d, x)]
  d.deltaval = [d.deltaval x]


  for i = 1:d.N+1
    grad_f[i] = (d.xpts[:,d.tauIdx[i+1]]'*d.C[:, :, i]*d.xpts[:,d.tauIdx[i+1]])[1]
  end

end






function MathProgBase.eval_hesslag(d::linSTOev, H, x, sigma, mu )

  # Increase counter of Hessian evaluations
  d.nhesseval += 1

  # Check if the matrices have already been precomputed
  if d.deltahess_prev!= x
    precompMatHessJ!(d, x)
    d.deltahess_prev[:] = x
  end

  Htemp = zeros(d.N+1, d.N+1)


  for i = 1:d.N+1
    Htemp[i, i] = 2*(d.xpts[:, d.tauIdx[i+1]]'*d.C[:, :, i]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
  end

  for j = 2:d.N+1
      for i = 1:j-1
      Htemp[j, i] = 2*(d.xpts[:, d.tauIdx[j+1]]'*d.C[:, :, j]*d.Phi[:, :, i+1, j+1]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
    end
  end

  Htemp *= sigma

  # Assign Elements of Hessian
  H[:] = Htemp[d.IndTril]
end



function MathProgBase.eval_hesslag(d::nlinSTOev, H, x, sigma, mu )

  # Increase counter of Hessian evaluations
  d.nhesseval += 1

  # Check if the matrices have already been precomputed
  if d.deltahess_prev!= x
    precompMatHessJ!(d, x)
    d.deltahess_prev[:] = x
  end

  Htemp = zeros(d.N+1, d.N+1)

  for i = 1:d.N+1
    Htemp[i, i] = 2*(d.xpts[:, d.tauIdx[i+1]]'*d.C[:, :, i]*d.A[:, :, d.tauIdx[i+1]-1]*d.xpts[:, d.tauIdx[i+1]])[1]
  end

  for j = 2:d.N+1
      for i = 1:j-1
      Htemp[j, i] = 2*(d.xpts[:, d.tauIdx[j+1]]'*d.C[:, :, j]*d.Phi[:, :, i+1, j+1]*d.A[:, :, d.tauIdx[i+1]-1]*d.xpts[:, d.tauIdx[i+1]])[1]
    end
  end

  H[:] = sigma*Htemp[d.IndTril]
end





# Constraints
function MathProgBase.eval_g(d::STOev, g, x)
  g[:] = d.Ag*x
end


function MathProgBase.eval_jac_g(d::STOev, J, x)
  J[:] = d.Vg
end


function MathProgBase.jac_structure(d::STOev)
  return d.Ig, d.Jg
end


# Hessian structure for Linear and Nonlinear Systems
function MathProgBase.hesslag_structure(d::STOev)
  return d.Itril, d.Jtril
end


# Linearize Dynamics
function linearizeDyn(nonlin_dyn::Function, nonlin_dyn_deriv::Function, x::Array{Float64,1}, u::Array{Float64,1})

  f = nonlin_dyn(x,u)
  df = nonlin_dyn_deriv(x, u)

  A = [df f-df*x; zeros(1, length(x)+1)]

end


function mergeSortFindIndex(tgrid::Array{Float64, 1}, tau::Array{Float64,1})

  ngrid = length(tgrid)
  N = length(tau)
  tauIdx = Array(Int, N+2); tauIdx[1] = 1; tauIdx[end]= N + ngrid
  tgridIdx = Array(Int, ngrid); tgridIdx[1] = 1; tgridIdx[end]= N + ngrid


  # Create merged and sorted time vector with grid and switching times
  ttemp = vcat(tgrid, tau)  # Concatenate grid and tau vector
  tidxtemp = sortperm(ttemp)  # Find permutation vector to sort ttemp
  tvec = ttemp[tidxtemp]    # Create full sorted tvec


  # # Create index of the tau vector elements inside tvec
  for i = 1:N
    tauIdx[i+1] = findfirst(tidxtemp, ngrid + i)
  end

  # # Create index of the tgrid vector elements inside tvec
  for i = 1:ngrid
    tgridIdx[i] = findfirst(tidxtemp, i)
  end

  return tvec, tauIdx, tgridIdx

end


# Propagate dynamic for linear STO
function propagateDynamics!(d::linSTOev, x::Array{Float64,1})

  # Get positive delta
  x = max(x, 0)

  # Get switching times from delta (tfdelta is the final time we get from the current delta vector)
  tau, d.tfdelta = delta2tau(x, d.t0)


  # Create grid from t0 to tfdelta
  d.tgrid[end] = d.tfdelta
  d.tgrid = min(d.tgrid, d.tgrid[end])


  # Create merged and sorted time vector with grid and switching times
  d.tvec, d.tauIdx, d.tgridIdx = mergeSortFindIndex(d.tgrid, tau)


  # Get complete delta vector with all intervals
  d.deltacomplete = tau2delta(d.tvec[2:end-1], d.t0, d.tfdelta)

  # The derivative checker in IPOPT will fail in computing the numerical derivatives because of the numerical issues in going to tau formulation and then back to delta formulation. To double check, please uncomment the following line in the case of only 2 points in the grid. The derivative checker should have no errors.

  # d.deltacomplete = x


  # Compute Exponentials required for the computations (for the complete grid)
  Aidx = 1  # Initialize index for current A mode

  # Compute Matrix Exponentials
  for i = 1:d.N+d.ngrid-1  # Iterate over all grid (incl sw times)


    # Verify which mode input applies
    if Aidx <= d.N
      if i>= d.tauIdx[Aidx + 1]
        Aidx += 1
      end
    end

    # Get temporary Matrix relative to current Aidx
    if d.isDiag[Aidx]  # Diagonalizable Matrix -> Compute Fast Exponential
      tempMat = real(d.V[:, :, Aidx]*diagm(exp(d.D[:,Aidx]*d.deltacomplete[i]))*d.invV[:, :, Aidx])
    else  # Nondiagonalizable Matrix -> Compute Standard Exponential
      # Compute Matrix Exponential
      tempMat = expm([-d.A[:, :, i]'  d.Q;
                      zeros(d.nx, d.nx)  d.A[:, :, i]]*d.deltacomplete[i])
    end



    # Assign \mathcal{E}_i
    d.expMat[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]

    # Assign M_k
    d.M[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]'*tempMat[1:d.nx, d.nx+1:2*d.nx]

    # Compute next state
    d.xpts[:, i+1] = d.expMat[:,:,i]*d.xpts[:, i]


  end

end




# Propagate dynamics for nonlinear STO
function propagateDynamics!(d::nlinSTOev, x::Array{Float64,1})


  # Get positive delta
  x = max(x, 0)

  # Get switching times from delta
  tau, d.tfdelta = delta2tau(x, d.t0)


  # Create grid from t0 to tfdelta
  d.tgrid[end] = d.tfdelta
  d.tgrid = min(d.tgrid, d.tgrid[end])


  # Fit switching times withing the grid
  d.tvec, d.tauIdx, d.tgridIdx = mergeSortFindIndex(d.tgrid, tau)


  # Get complete delta vector with all intervals
  d.deltacomplete = tau2delta(d.tvec[2:end-1], d.t0, d.tfdelta)

  # Propagate Dynamics and Compute Exponentials over the whole grid
  uIdx = 1  # Initialize index for current u

  # Compute Matrix Exponentials
  for i =1:d.N+d.ngrid-1  # Iterate over all grid (incl sw times)

    # Verify which U input applies
    if uIdx <= d.N
      if i>= d.tauIdx[uIdx + 1]
        uIdx += 1
      end
    end


    # Linearize Dynamics
    d.A[:,:,i] = linearizeDyn(d.nonlin_dyn, d.nonlin_dyn_deriv, d.xpts[1:end-1,i], d.uvec[:, uIdx])

    # Compute Matrix Exponential
    tempMat = expm([-d.A[:, :, i]'  d.Q;
                    zeros(d.nx, d.nx)  d.A[:, :, i]]*d.deltacomplete[i])

    # Assign \mathcal{E}_i
    d.expMat[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]

    # Assign M_k
    d.M[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]'*tempMat[1:d.nx, d.nx+1:2*d.nx]

    # Compute state at the next grid point
    d.xpts[:, i+1] = d.expMat[:, :, i]*d.xpts[:, i]

  end




end

"Low level function for computing S matrices"
function computeSlw!(d::STOev)
  d.S[:,:,d.N+d.ngrid] = d.E
  for i = d.N+d.ngrid-1:-1:1
    d.S[:,:,i] = d.M[:,:,i] + d.expMat[:,:,i]'*d.S[:,:,i+1]*d.expMat[:,:,i]
  end
end


"Low level function for computing C matrices for linear dynamics"
function computeClw!(d::linSTOev)
  for i = 1:d.N+1
    d.C[:, :, i] = d.Q + d.A[:, :, i]'*d.S[:, :, d.tauIdx[i+1]] + d.S[:, :, d.tauIdx[i+1]]*d.A[:, :, i]
  end
end

"Low level function for computing C matrices for nonlinear dynamics"
function computeClw!(d::nlinSTOev)
  for i = 1:d.N+1
    d.C[:, :, i] = d.Q + d.A[:, :, d.tauIdx[i+1]-1]'*d.S[:, :, d.tauIdx[i+1]] + d.S[:, :, d.tauIdx[i+1]]*d.A[:, :, d.tauIdx[i+1]-1]
  end
end



"Low level function for computing Phi matrices"
function computePhilw!(d::STOev)
  for i = 1 : d.N + 2
    d.Phi[:, :, i, i] = eye(d.nx)  # Identity Matrix to Start
    for l = i:d.N + 1  # Iterate over all successive Phi matrices
      d.Phi[:, :, i, l+1] = d.Phi[:, :, i, l]  # Initialize with previous matrix
      for j=d.tauIdx[l]:d.tauIdx[l+1]-1  # Iterate over all points between switching times
        d.Phi[:, :, i, l+1] = d.expMat[:, :, j] * d.Phi[:, :, i, l+1]
      end
    end
  end
end
