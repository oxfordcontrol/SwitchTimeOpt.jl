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


#  Function for Nonlinear STO
function precompMatrices!(d::nlinSTOev, x)
  x = [d.t0; x; d.tf]  # Create Vector with initial and final times

  # Compute Matrix Exponentials
  for i = 1:d.N+1

    # Linearize Dynamics
    d.A[:,:,i] = linearizeDyn(d.nonlin_dyn, d.nonlin_dyn_deriv,       d.xpts[1:end-1,i], d.uvec[:,i])

    # Compute Matrix Exponential
    tempMat = expm([-d.A[:, :, i]'  d.Q;
                    zeros(d.nx, d.nx)  d.A[:, :, i]]*(x[i+1] - x[i]))
    # tempMat = real(d.V[:, :, i]*diagm(exp(d.D[:,i]*(x[i+1] - x[i])))*d.invV[:, :, i])

    # Assign \mathcal{E}_i
    d.expMat[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]

    # Assign M_k
    d.M[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]'*tempMat[1:d.nx, d.nx+1:2*d.nx]

  end


  # Compute State Transition Matrices
  for i = 1:d.N+1
    d.Phi[:,:,i,i] = eye(d.nx)  # Identity Matrix to Start
    for j = i+1:d.N+1
      d.Phi[:,:,i,j] = d.expMat[:,:,j-1]*d.Phi[:,:,i, j-1]
    end
  end

  # Compute States at Switching Times
  d.xpts[:, 1] = d.x0
  for i = 1:d.N
    d.xpts[:, i+1] = d.expMat[:,:,i]*d.xpts[:, i]
  end

  # Compute P integrals
  d.P[:,:,d.N+1] = d.M[:,:, d.N+1]
  for i = d.N:-1:1
    d.P[:,:,i] = d.M[:,:,i] + d.expMat[:,:,i]'*d.P[:,:,i+1]*d.expMat[:,:,i]
  end

end


function precompMatrices!(d::linSTOev, x)
  x = [d.t0; x; d.tf]  # Create Vector with initial and final times

  # Compute Matrix Exponentials
  for i = 1:d.N+1

    # Compute Matrix Exponential to
    # tempMat = expm([-d.A[:, :, i]'  d.Q;
                    # zeros(d.nx, d.nx)  d.A[:, :, i]]*(x[i+1] - x[i]))
    tempMat = real(d.V[:, :, i]*diagm(exp(d.D[:,i]*(x[i+1] - x[i])))*d.invV[:, :, i])

    # Assign \mathcal{E}_i
    d.expMat[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]

    # Assign M_k
    d.M[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]'*tempMat[1:d.nx, d.nx+1:2*d.nx]

  end


  # Compute State Transition Matrices
  for i = 1:d.N+1
    d.Phi[:,:,i,i] = eye(d.nx)  # Identity Matrix to Start
    for j = i+1:d.N+1
      d.Phi[:,:,i,j] = d.expMat[:,:,j-1]*d.Phi[:,:,i, j-1]
    end
  end

  # Compute States at Switching Times
  d.xpts[:, 1] = d.x0
  for i = 1:d.N
    d.xpts[:, i+1] = d.expMat[:,:,i]*d.xpts[:, i]
  end

  # Compute P integrals
  d.P[:,:,d.N+1] = d.M[:,:, d.N+1]
  for i = d.N:-1:1
    d.P[:,:,i] = d.M[:,:,i] + d.expMat[:,:,i]'*d.P[:,:,i+1]*d.expMat[:,:,i]
  end

end




# Evaluate Function
function MathProgBase.eval_f(d::STOev, x)
  # Check if the matrices have already been precomputed
  if d.prev_tau != x
      precompMatrices!(d, x) # Precompute Matrices and store them in d
      d.prev_tau[:] = x  # Update current tau
  end
  J = (0.5*d.x0'*d.P[:,:,1]*d.x0)[1]

end



function MathProgBase.eval_grad_f(d::STOev, grad_f, x)
  # Check if the matrices have already been precomputed
  if d.prev_tau != x
      precompMatrices!(d, x) # Precompute Matrices and store them in d
      d.prev_tau[:] = x  # Update current tau
  end

  for i = 2:d.N+1
    grad_f[i-1] = (d.xpts[:,i]'*d.P[:,:,i]*(d.A[:,:,i-1] - d.A[:,:,i])*d.xpts[:,i])[1]
  end

end

function MathProgBase.eval_hesslag(d::STOev, H, x, sigma, mu )

  # Check if the matrices have already been precomputed
  if d.prev_tau != x
      precompMatrices!(d, x) # Precompute Matrices and store them in d
      d.prev_tau[:] = x  # Update current tau
  end

  Htemp = zeros(d.N, d.N)

  # Evaluate double derivative with respect to the same switching instants
  for i = 2:d.N+1
    Htemp[i-1, i-1] = (d.xpts[:, i]'*(-d.Q*(d.A[:, :, i-1] - d.A[:, :, i]) + (d.A[:, :, i-1] - d.A[:, :, i])'*d.P[:, :, i]*(d.A[:, :, i-1] - d.A[:, :, i]) + d.P[:, :, i]*(d.A[:, :, i-1]^2 + d.A[:, :, i]^2 - 2*d.A[:, :, i]*d.A[:,:,i-1]))*d.xpts[:, i])[1]
  end

  # Evaluate double derivative for mixed switching instants
  for i = 2:d.N
    for j = i+1:d.N+1
      Htemp[j-1, i-1] = (d.xpts[:, i]'*(d.A[:, :, i-1] - d.A[:, :, i])'*d.Phi[:, :, i, j]'*(d.P[:, :, j]*(d.A[:, :, j-1] - d.A[:, :, j]) + (d.A[:, :, j-1] - d.A[:, :, j])'*d.P[:, :, j])*d.xpts[:, j])[1]
    end
  end

  # H[:,:] = sigma*H
  # _, _, H[:] = findnz(sigma*Htemp)
  # Ind = find(tril(ones(d.N, d.N)))
  # H[:] = sigma*Htemp[Ind]
  H[:] = sigma*Htemp[d.IndTril]
end


function MathProgBase.eval_g(d::STOev, g, x)
  # Mg = [-[eye(d.N-1) zeros(d.N-1)] + [zeros(d.N-1) eye(d.N-1)]; 1 zeros(1,d.N-1); zeros(1,d.N-1) -1]
  #
  # g[:] = Mg*x + [zeros(d.N-1); -d.t0; d.tf]

  g[:] = d.Ag*x - d.bg

end


function MathProgBase.eval_jac_g(d::STOev, J, x)
  # Mg = [-[eye(d.N-1) zeros(d.N-1,1)] + [zeros(d.N-1,1) eye(d.N-1)]; 1 zeros(1,d.N-1); zeros(1,d.N-1) -1]
  # _, _, V = findnz(Mg)
  # J[:] = V
  J[:] = d.Vg
end


function MathProgBase.jac_structure(d::STOev)
  #   Mg = sparse([-[eye(d.N-1) zeros(d.N-1)] + [zeros(d.N-1) eye(d.N-1)]; 1 zeros(1,d.N-1); zeros(1,d.N-1) -1])
  #   Iind, Jind, _ = findnz(Mg)
  # return Iind, Jind
  return d.Ig, d.Jg
end

function MathProgBase.hesslag_structure(d::STOev)
  # # Define upper triangular Matrix
  # Iind, Jind, _ = findnz(tril(ones(d.N, d.N)))
  # return Iind, Jind
  return d.Itril, d.Jtril
end


# Linearize Dynamics
function linearizeDyn(nonlin_dyn::Function, nonlin_dyn_deriv::Function, x::Array{Float64,1}, u::Array{Float64,1})

  f = nonlin_dyn(x,u)
  df = nonlin_dyn_deriv(x, u)

  A = [df f-df*x; zeros(1, length(x)+1)]

end
