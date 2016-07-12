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


#  Function for Linear STO
function precompMatrices!(d::linSTOev, x)
  # x = [d.t0; x; d.tf]  # Create Vector with initial and final times

  #-------------------------------------------------------
  # Construct grid with new delta vector x
  #-------------------------------------------------------

  # Get switching times from delta
  tau = delta2tau(x, d.t0, d.tf)

  # Propagate dynamic for new switching times
  propagateDynamics!(d, tau)


  # The derivative checker in IPOPT will fail in computing the numerical derivatives because of the numerical issues in going to tau formulation and then back to delta formulation. To double check, please uncomment the following line in the case of only 2 points in the grid. The derivative checker should have no errors.
  # d.deltacomplete = x

  # # Create merged and sorted time vector with grid and switching times
  # ttemp = vcat(d.tgrid, tau)  # Concatenate grid and tau vector
  # tidxtemp = sortperm(ttemp)  # Find permutation vector to sort ttemp
  # d.tvec = ttemp[tidxtemp]    # Create full sorted tvec
  #
  # # # Create index of the tau vector elements inside tvec
  # for i = 1:d.N
  #   d.tauIdx[i+1] = findfirst(tidxtemp, d.ngrid + i)
  # end
  #
  #
  # # Create merged and sorted time vector with grid and switching times
  # d.tvec, d.tauIdx = mergeSortFindIndex(d.tgrid, tau)
  # # tvectest, tauIdxtest = mergeSortFindIndex(d.tgrid, tau)
  #
  #
  #
  # # # Create merged and sorted time vector with grid and switching times
  # # d.tvec = sort(vcat(d.tgrid, tau))
  # #
  # # # Create index of the tau vector elements inside tvec
  # # for i = 1:d.N
  # #   d.tauIdx[i+1] = findfirst(d.tvec, tau[i])  # i+1 because tau0ws
  # #
  # #   # # Check if tau vector is at the end (IS IT NECESSARY?)
  # #   # if d.tauIdx[i+1] == d.N + d.ngrid
  # #   #   d.tauIdx[i+1] -= 1
  # #   # end
  # #
  # # end
  # #
  # # if d.tauIdx != tauIdxtest
  # #   @printf("\nError in two functions!\n\n")
  # #   @printf("tau = "); show(tau); @printf("\n")
  # #   @printf("tgrid = "); show(d.tgrid); @printf("\n")
  # #   @printf("tvec = "); show(d.tvec); @printf("\n")
  # #   @printf("tvectest = "); show(tvectest); @printf("\n")
  # #   @printf("tauIdx = "); show(d.tauIdx); @printf("\n")
  # #   @printf("tauIdxtest = "); show(tauIdxtest); @printf("\n")
  # # end
  #
  #
  # # Get complete delta vector with all intervals
  # d.deltacomplete = tau2delta(d.tvec[2:end-1], d.t0, d.tf)
  #
  # # The derivative checker in IPOPT will fail in computing the numerical derivatives because of the numerical issues in going to tau formulation and then back to delta formulation. To double check, please uncomment the following line in the case of only 2 points in the grid. The derivative checker should have no errors.
  #
  #
  # # d.deltacomplete = x
  #
  #
  # #-----------------------------------------------------------------------------
  # # Compute Exponentials required for the computations (for the complete grid)
  # #-----------------------------------------------------------------------------
  #
  # Aidx = 1  # Initialize index for current A mode
  #
  # # Compute Matrix Exponentials
  # for i = 1:d.N+d.ngrid-1  # Iterate over all grid (incl sw times)
  #
  #
  #   # Verify which mode input applies
  #   if Aidx <= d.N
  #     if i>= d.tauIdx[Aidx + 1]
  #       Aidx += 1
  #     end
  #   end
  #
  #   # Compute Matrix Exponential to
  #   # tempMat = expm([-d.A[:, :, i]'  d.Q;
  #                   # zeros(d.nx, d.nx)  d.A[:, :, i]]*x[i])
  #   # tempMat = real(d.V[:, :, i]*diagm(exp(d.D[:,i]*x[i]))*d.invV[:, :, i])
  #
  #   # # Get temporary Matrix relative to current Aidx
  #   if d.isDiag[Aidx]  # Diagonalizable Matrix -> Compute Fast Exponential
  #     tempMat = real(d.V[:, :, Aidx]*diagm(exp(d.D[:,Aidx]*d.deltacomplete[i]))*d.invV[:, :, Aidx])
  #   else  # Nondiagonalizable Matrix -> Compute Standard Exponential
  #     # Compute Matrix Exponential
  #     tempMat = expm([-d.A[:, :, i]'  d.Q;
  #                     zeros(d.nx, d.nx)  d.A[:, :, i]]*d.deltacomplete[i])
  #   end
  #
  #
  #
  #   # Assign \mathcal{E}_i
  #   d.expMat[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]
  #
  #   # Assign M_k
  #   d.M[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]'*tempMat[1:d.nx, d.nx+1:2*d.nx]
  #
  # end



  #------------------------------------------------------------
  # Compute State Transition Matrices (for the complete grid)
  #------------------------------------------------------------
  for i = 1:d.N+d.ngrid
    d.Phi[:,:,i,i] = eye(d.nx)  # Identity Matrix to Start
    for j = i+1:d.N+d.ngrid
      d.Phi[:,:,i,j] = d.expMat[:,:,j-1]*d.Phi[:,:,i, j-1]
    end
  end

  # # Compute States at Switching Times
  # d.xpts[:, 1] = d.x0
  # for i = 1:d.N+d.ngrid-1
  #   d.xpts[:, i+1] = d.expMat[:,:,i]*d.xpts[:, i]
  # end


  #-----------------------------------------------------------------------
  # Compute matrices S and C for gradient and hessian
  #-----------------------------------------------------------------------

  # Compute S matrices for the whole grid
  d.S[:,:,d.N+d.ngrid] = d.Qf
  # d.S[:,:,d.N+d.ngrid-1] = d.M[:,:, end] + d.expMat[:,:, end]'*d.Qf*d.expMat[:,:, end]

  for i = d.N+d.ngrid-1:-1:1
    d.S[:,:,i] = d.M[:,:,i] + d.expMat[:,:,i]'*d.S[:,:,i+1]*d.expMat[:,:,i]
  end

  # Compute C Matrices for every switching time
  # The subsequent S[i+1]*A[i] means the indexing of the complete grid at the interval i
  #

  for i = 1:d.N+1
    d.C[:, :, i] = d.Q + d.A[:, :, i]'*d.S[:, :, d.tauIdx[i+1]] + d.S[:, :, d.tauIdx[i+1]]*d.A[:, :, i]
  end


  #-----------------------------------------------------------------------
  # Compute Constraints Jacobian
  #-----------------------------------------------------------------------
  # if d.ncons!=0
  #
  #   # Compute Jacobian
  #   Jac_temp = zeros(1+d.ncons, d.N+1)
  #   Jac_temp[1,:] = d.gsum  # Constraint on the sum of Variables
  #
  #   # Construct jacobian (Constraint on last stage)
  #   for l = 1:d.ncons # Iterate over Constraints
  #     for i = 1:d.N+1 # Iterate over Variables
  #       Jac_temp[1+l, i] = (d.Ac[l, :]*d.Phi[:, :, d.tauIdx[i+1], end]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
  #     end
  #   end
  #
  #   d.Vg = Jac_temp[:]
  #   # Easy to check where index k of the state constraint is before or after the current tau. Just check whether k is greater or lower than tauIdx. (TODO!)
  #
  #
  #
  # end








  #
  #
  # # Compute State Transition Matrices
  # for i = 1:d.N+2
  #   d.Phi[:,:,i,i] = eye(d.nx)  # Identity Matrix to Start
  #   for j = i+1:d.N+2
  #     d.Phi[:,:,i,j] = d.expMat[:,:,j-1]*d.Phi[:,:, i, j-1]
  #   end
  # end
  #
  # # Compute States at Switching Times
  # d.xpts[:, 1] = d.x0
  # for i = 1:d.N+1
  #   d.xpts[:, i+1] = d.expMat[:,:,i]*d.xpts[:, i]
  # end
  #
  # # Compute S matrices
  # d.S[:,:,d.N+1] = d.M[:,:, d.N+1] + d.expMat[:,:,d.N+1]'*d.Qf*d.expMat[:,:,d.N+1]
  #
  # for i = d.N:-1:1
  #   d.S[:,:,i] = d.M[:,:,i] + d.expMat[:,:,i]'*d.S[:,:,i+1]*d.expMat[:,:,i]
  # end
  #
  # # Compute C matrices
  # for i = 1:d.N
  #   d.C[:, :, i] = d.Q + d.A[:,:,i]'*d.S[:,:,i+1] + d.S[:,:,i+1]*d.A[:, :, i]
  # end
  # d.C[:, :, d.N+1] = d.Q + d.A[:,:, d.N+1]'*d.Qf + d.Qf*d.A[:,:, d.N+1] # S[:, :, N+2] = Q_f
end






#  Function for Nonlinear STO
function precompMatrices!(d::nlinSTOev, x)

  #-------------------------------------------------------
  # Construct grid with new delta vector x
  #-------------------------------------------------------

  # Get switching times from delta
  tau = delta2tau(x, d.t0, d.tf)

  # Propagate Dynamics to compute matrix exponentials and states at the switching times
  propagateDynamics!(d, tau)

  # d.tvec, d.tauIdx = mergeSortFindIndex(d.tgrid, tau)
  #
  # # d.tvec = sort(vcat(d.tgrid, tau))
  # #
  # # # Create index of the tau vector elements inside tvec
  # # for i = 1:d.N
  # #   # d.tauIdx[i+1] = findfirst(d.tvec, tau[i])  # i+1 because tau0ws
  # #   d.tauIdx[i+1] = findnext(d.tvec, tau[i], d.tauIdx[i])  # find next (greater than last tauIdx)
  # #
  # #   # # Check if tau vector is at the end (IS IT NECESSARY?)
  # #   # if d.tauIdx[i+1] == d.N + d.ngrid
  # #   #   d.tauIdx[i+1] -= 1
  # #   # end
  # #
  # # end
  #
  #
  # # Get complete delta vector with all intervals
  # d.deltacomplete = tau2delta(d.tvec[2:end-1], d.t0, d.tf)
  #
  # #-----------------------------------------------------------------------------
  # # Compute Exponentials required for the computations (for the complete grid)
  # #-----------------------------------------------------------------------------
  #
  # uIdx = 1  # Initialize index for current u
  #
  # # Compute Matrix Exponentials
  # for i =1:d.N+d.ngrid-1  # Iterate over all grid (incl sw times)
  #
  #   # Verify which U input applies
  #   if uIdx <= d.N
  #     if i>= d.tauIdx[uIdx + 1]
  #       uIdx += 1
  #     end
  #   end
  #
  #
  #   # Linearize Dynamics
  #   d.A[:,:,i] = linearizeDyn(d.nonlin_dyn, d.nonlin_dyn_deriv, d.xpts[1:end-1,i], d.uvec[:,uIdx])
  #   # d.A[:,:,i] = linearizeDyn(d.nonlin_dyn, d.nonlin_dyn_deriv, d.xpts[1:end-1,i+1], d.uvec[:,uIdx])  # try linearization at next step
  #
  #   # Compute Matrix Exponential
  #   tempMat = expm([-d.A[:, :, i]'  d.Q;
  #                   zeros(d.nx, d.nx)  d.A[:, :, i]]*d.deltacomplete[i])
  #
  #
  #   # Assign \mathcal{E}_i
  #   d.expMat[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]
  #
  #   # Assign M_k
  #   d.M[:,:,i] = tempMat[d.nx+1:2*d.nx, d.nx+1:2*d.nx]'*tempMat[1:d.nx, d.nx+1:2*d.nx]
  #
  # end

  #------------------------------------------------------------
  # Compute State Transition Matrices (for the complete grid)
  #------------------------------------------------------------
  for i = 1:d.N+d.ngrid
    d.Phi[:,:,i,i] = eye(d.nx)  # Identity Matrix to Start
    for j = i+1:d.N+d.ngrid
      d.Phi[:,:,i,j] = d.expMat[:,:,j-1]*d.Phi[:,:,i, j-1]
    end
  end

  # # Compute States at Switching Times
  # d.xpts[:, 1] = d.x0
  # for i = 1:d.N+d.ngrid-1
  #   d.xpts[:, i+1] = d.expMat[:,:,i]*d.xpts[:, i]
  # end


  #-----------------------------------------------------------------------
  # Compute matrices S and C for gradient and hessian
  #-----------------------------------------------------------------------

  # Compute S matrices for the whole grid
  d.S[:,:,d.N+d.ngrid] = d.Qf
  # d.S[:,:,d.N+d.ngrid-1] = d.M[:,:, end] + d.expMat[:,:, end]'*d.Qf*d.expMat[:,:, end]

  for i = d.N+d.ngrid-1:-1:1
    d.S[:,:,i] = d.M[:,:,i] + d.expMat[:,:,i]'*d.S[:,:,i+1]*d.expMat[:,:,i]
  end


  # # Compute S matrices
  # Stemp = d.M[:,:, end] + d.expMat[:, :, end]'*d.Qf*d.expMat[:, :, end]
  # for i = d.N+d.ngrid-2:-1:1
  #
  #   # Update Temporary Matrix
  #   Stemp = d.M[:,:,i] + d.expMat[:,:,i]'*Stemp*d.expMat[:,:,i]
  #
  #   # Check i index corresponds to one of the switching times
  #   if (tauSwIdx = findfirst(d.tauIdx, i))!=0
  #     d.S[:, :, tauSwIdx+1] = Stemp  # Add S matrix corresponding to the SW time (+1 because S has one more matrix than the number of switches (starts at time t0))
  #   end
  # end
  #
  # d.S[:, :, 1] = Stemp  # Add initial Matrix S


  # # Compute S matrices
  # d.S[:,:,d.N+1] = d.M[:,:, d.N+1] + d.expMat[:,:,d.N+1]'*d.Qf*d.expMat[:,:,d.N+1]
  #
  # for i = d.N:-1:1
  #   d.S[:,:,i] = d.M[:,:,i] + d.expMat[:,:,i]'*d.S[:,:,i+1]*d.expMat[:,:,i]
  # end



  # # Compute C matrices
  # d.C[:, :, 1] = d.Q + d.A[:,:,1]'*d.S[:,:,2] + d.S[:,:,2]*d.A[:, :, 1]
  #
  # for i = 2:d.N
  #   d.C[:, :, i] = d.Q + d.A[:,:,d.tauIdx[i-1]]'*d.S[:,:,i+1] + d.S[:,:,i+1]*d.A[:, :, d.tauIdx[i-1]]
  # end
  # d.C[:, :, d.N+1] = d.Q + d.A[:,:, d.tauIdx[d.N]]'*d.Qf + d.Qf*d.A[:,:, d.tauIdx[d.N]] # S[:, :, N+ngrid] = Q_f


  # Compute C Matrices for every switching time
  # The subsequent S[i+1]*A[i] means the indexing of the complete grid at the interval i
  #

  for i = 1:d.N+1
    d.C[:, :, i] = d.Q + d.A[:, :, d.tauIdx[i+1]-1]'*d.S[:, :, d.tauIdx[i+1]] + d.S[:, :, d.tauIdx[i+1]]*d.A[:, :, d.tauIdx[i+1]-1]
  end

  # Backup working one
  # for i = 1:d.N+1
  #   d.C[:, :, i] = d.Q + d.A[:, :, d.tauIdx[i]]'*d.S[:, :, d.tauIdx[i]+1] + d.S[:, :, d.tauIdx[i]+1]*d.A[:, :, d.tauIdx[i]]
  # end








  # d.C[:, :, 1] = d.Q + d.A[:, :, 1]'*d.S[:, :, 2] + d.S[:, :, 2]*d.A[:, :, 1]
  #
  # for i = 2:d.N+1
  #   d.C[:, :, i] = d.Q + d.A[:, :, d.tauIdx[i-1]]'*d.S[:, :, d.tauIdx[i-1]+1] + d.S[:, :, d.tauIdx[i-1]+1]*d.A[:, :, d.tauIdx[i-1]]
  # end


  # # Compute C matrices
  # for i = 1:d.N
  #   d.C[:, :, i] = d.Q + d.A[:,:,i]'*d.S[:,:,i+1] + d.S[:,:,i+1]*d.A[:, :, i]
  # end
  # d.C[:, :, d.N+1] = d.Q + d.A[:,:, d.N+1]'*d.Qf + d.Qf*d.A[:,:, d.N+1] # S[:, :, N+2] = Q_f

end




# Evaluate Cost Function
function MathProgBase.eval_f(d::STOev, x)
  # Check if the matrices have already been precomputed
  if d.prev_delta != x
      precompMatrices!(d, x) # Precompute Matrices and store them in d
      d.prev_delta[:] = x  # Update current tau
  end
  J = 0.5*(d.x0'*d.S[:,:,1]*d.x0)[1]
end



# # Evaluate Gradient for Nonlinear Dynamics
# function MathProgBase.eval_grad_f(d::nlinSTOev, grad_f, x)
#   # Check if the matrices have already been precomputed
#   if d.prev_delta != x
#       precompMatrices!(d, x) # Precompute Matrices and store them in d
#       d.prev_delta[:] = x  # Update current tau
#   end
#
#   # for i = 2:d.N+1
#   #   grad_f[i-1] = (d.xpts[:,i]'*d.P[:,:,i]*(d.A[:,:,i-1] - d.A[:,:,i])*d.xpts[:,i])[1]
#   # end
#
#   # Working without grid
#   # for i = 1:d.N+1
#   #   grad_f[i] = 0.5*(d.xpts[:,i+1]'*d.C[:, :, i]*d.xpts[:,i+1])[1]
#   # end
#
#
#   for i = 1:d.N+1
#     grad_f[i] = 0.5*(d.xpts[:,d.tauIdx[i+1]]'*d.C[:, :, i]*d.xpts[:,d.tauIdx[i+1]])[1]
#   end
#
#   # grad_f[1] = 0.5*(d.xpts[:,2]'*d.C[:, :, 1]*d.xpts[:,2])[1]
#   # for i = 2:d.N+1
#   #   grad_f[i] = 0.5*(d.xpts[:,d.tauIdx[i-1]+1]'*d.C[:, :, i]*d.xpts[:,d.tauIdx[i-1]+1])[1]
#   # end
#
# end



# Evaluate Gradient
function MathProgBase.eval_grad_f(d::STOev, grad_f, x)
  # Check if the matrices have already been precomputed
  if d.prev_delta != x
      precompMatrices!(d, x) # Precompute Matrices and store them in d
      d.prev_delta[:] = x  # Update current tau
  end

  # for i = 2:d.N+1
  #   grad_f[i-1] = (d.xpts[:,i]'*d.P[:,:,i]*(d.A[:,:,i-1] - d.A[:,:,i])*d.xpts[:,i])[1]
  # end

  # Working without grid
  # for i = 1:d.N+1
  #   grad_f[i] = 0.5*(d.xpts[:,i+1]'*d.C[:, :, i]*d.xpts[:,i+1])[1]
  # end


  for i = 1:d.N+1
    grad_f[i] = 0.5*(d.xpts[:,d.tauIdx[i+1]]'*d.C[:, :, i]*d.xpts[:,d.tauIdx[i+1]])[1]
  end

  # grad_f[1] = 0.5*(d.xpts[:,2]'*d.C[:, :, 1]*d.xpts[:,2])[1]
  # for i = 2:d.N+1
  #   grad_f[i] = 0.5*(d.xpts[:,d.tauIdx[i-1]+1]'*d.C[:, :, i]*d.xpts[:,d.tauIdx[i-1]+1])[1]
  # end

end






function MathProgBase.eval_hesslag(d::linSTOev, H, x, sigma, mu )

  # Check if the matrices have already been precomputed
  if d.prev_delta != x
      precompMatrices!(d, x) # Precompute Matrices and store them in d
      d.prev_delta[:] = x  # Update current tau
  end

  Htemp = zeros(d.N+1, d.N+1)


  for i = 1:d.N+1
    Htemp[i, i] = (d.xpts[:, d.tauIdx[i+1]]'*d.C[:, :, i]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
  end

  for j = 2:d.N+1
    # for j = i+1:d.N+1
      for i = 1:j-1
      Htemp[j, i] = (d.xpts[:, d.tauIdx[j+1]]'*d.C[:, :, j]*d.Phi[:, :, d.tauIdx[i+1], d.tauIdx[j+1]]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
    end
  end

  # H[:] = sigma*Htemp[d.IndTril]
  Htemp *= sigma

  # Add Stagewise Constraints
  if d.ncons!=0
    for k = 2:d.ngrid-1  # iterate over grid constraints
      for l = 1:d.ncons   # iterate over elements of the stagewise constraints
        # Initialize new empty Hgtemp
        Hgtemp = zeros(d.N+1, d.N+1)

        # double derivatives with respect to the same interval
        for i = 1:d.N+1
          if d.tauIdx[i+1] < d.tgridIdx[k]  # Change contraint only if not zero.
            Hgtemp[i, i] = (d.Ac[l, :]*d.Phi[:, :, d.tauIdx[i+1], d.tgridIdx[k]]*(d.A[:, :, i]^2)*d.xpts[:, d.tauIdx[i+1]])[1]
          end
        end

        # double derivatives with respect to different intervals
        for j = 2:d.N+1
          for i = 1:j-1
            if d.tauIdx[j+1] < d.tgridIdx[k]  # Change contraint only if not zero.
              Hgtemp[j, i] = (d.Ac[l, :]*d.Phi[:, :, d.tauIdx[j+1], d.tgridIdx[k]]*d.A[:, :, j]*d.Phi[:, :, d.tauIdx[i+1], d.tauIdx[j+1]]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
            end
          end
        end


        Htemp += mu[1 + (k-2)*d.ncons+l]*Hgtemp  # Add 1 for linear sum constraint

      end  # End l for
    end  # End k for
  end

  # Final Stage Constraints
  if d.nconsf!=0 # If there are any constraints in the problem
    # Add hessian constraints on states
    Hgtemp = zeros(d.N+1, d.N+1)


    for l = 1:d.nconsf # Iterate over Constraints
      # Work on Hgtempf (hessian constraint l)

      for i = 1:d.N+1
        Hgtemp[i, i] = (d.Acf[l, :]*d.Phi[:, :, d.tauIdx[i+1], end]*(d.A[:, :, i]^2)*d.xpts[:, d.tauIdx[i+1]])[1]
      end

      for j = 2:d.N+1
          for i = 1:j-1
          Hgtemp[j, i] = (d.Acf[l, :]*d.Phi[:, :, d.tauIdx[j+1], end]*d.A[:, :, j]*d.Phi[:, :, d.tauIdx[i+1], d.tauIdx[j+1]]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
        end
      end

      Htemp += mu[1+(d.ngrid - 2)*d.ncons+l]*Hgtemp  # Add 1 for linear sum constraint
    end

  end


  # Assign Elements of Hessian
  H[:] = Htemp[d.IndTril]
end



 # # Old Hessian Stuff
 #  if d.ncons!=0 # If there are any constraints in the problem
 #    # Add hessian constraints on states
 #    Hgtempl = zeros(d.N+1, d.N+1)
 #
 #
 #    for l = 1:d.ncons # Iterate over Constraints
 #      # Work on Hgtempl (hessian constraint l)
 #
 #      for i = 1:d.N+1
 #        Hgtempl[i, i] = (d.Ac[l, :]*d.Phi[:, :, d.tauIdx[i+1], end]*(d.A[:, :, i]^2)*d.xpts[:, d.tauIdx[i+1]])[1]
 #      end
 #
 #      for j = 2:d.N+1
 #          for i = 1:j-1
 #          Hgtempl[j, i] = (d.Ac[l, :]*d.Phi[:, :, d.tauIdx[j+1], end]*d.A[:, :, j]*d.Phi[:, :, d.tauIdx[i+1], d.tauIdx[j+1]]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
 #        end
 #      end
 #
 #      # @printf("This look is active!")
 #
 #      Htemp += mu[1+l]*Hgtempl  # Add 1 for linear sum constraint
 #    end
 #
 #  end
  # Assign Elements of Hessian
   #  H[:] = Htemp[d.IndTril]

# end


function MathProgBase.eval_hesslag(d::nlinSTOev, H, x, sigma, mu )

  # Check if the matrices have already been precomputed
  if d.prev_delta != x
      precompMatrices!(d, x) # Precompute Matrices and store them in d
      d.prev_delta[:] = x  # Update current tau
  end

  Htemp = zeros(d.N+1, d.N+1)

  ## Evaluate double derivative with respect to the same switching instants
  # for i = 2:d.N+1
  #   Htemp[i-1, i-1] = (d.xpts[:, i]'*(-d.Q*(d.A[:, :, i-1] - d.A[:, :, i]) + (d.A[:, :, i-1] - d.A[:, :, i])'*d.P[:, :, i]*(d.A[:, :, i-1] - d.A[:, :, i]) + d.P[:, :, i]*(d.A[:, :, i-1]^2 + d.A[:, :, i]^2 - 2*d.A[:, :, i]*d.A[:,:,i-1]))*d.xpts[:, i])[1]
  # end
  #
  # # Evaluate double derivative for mixed switching instants
  # for i = 2:d.N
  #   for j = i+1:d.N+1
  #     Htemp[j-1, i-1] = (d.xpts[:, i]'*(d.A[:, :, i-1] - d.A[:, :, i])'*d.Phi[:, :, i, j]'*(d.P[:, :, j]*(d.A[:, :, j-1] - d.A[:, :, j]) + (d.A[:, :, j-1] - d.A[:, :, j])'*d.P[:, :, j])*d.xpts[:, j])[1]
  #   end
  # end

  # H[:,:] = sigma*H
  # _, _, H[:] = findnz(sigma*Htemp)
  # Ind = find(tril(ones(d.N, d.N)))
  # H[:] = sigma*Htemp[Ind]




  # Working with no grid
  # for i = 1:d.N+1
  #   Htemp[i, i] = (d.xpts[:, i+1]'*d.C[:, :, i]*d.A[:, :, i]*d.xpts[:, i+1])[1]
  # end
  #
  # for j = 2:d.N+1
  #   # for j = i+1:d.N+1
  #     for i = 1:j-1
  #     Htemp[j, i] = (d.xpts[:, j+1]'*d.C[:, :, j]*d.Phi[:, :, i+1, j+1]*d.A[:, :, i]*d.xpts[:, i+1])[1]
  #   end
  # end



  # Htemp[1, 1] = (d.xpts[:, 2]'*d.C[:, :, 1]*d.A[:, :, 1]*d.xpts[:, 2])[1]

  for i = 1:d.N+1
    Htemp[i, i] = (d.xpts[:, d.tauIdx[i+1]]'*d.C[:, :, i]*d.A[:, :, d.tauIdx[i+1]-1]*d.xpts[:, d.tauIdx[i+1]])[1]
    # Htemp[i, i] = 0.5*(d.xpts[:, d.tauIdx[i+1]]'*(d.A[:, :, d.tauIdx[i+1] - 1]'*d.C[:, :, i] + d.C[:, :, i]*d.A[:, :, d.tauIdx[i+1] - 1])*d.xpts[:, d.tauIdx[i+1]])[1]
  end

  for j = 2:d.N+1
    # for j = i+1:d.N+1
      for i = 1:j-1
      Htemp[j, i] = (d.xpts[:, d.tauIdx[j+1]]'*d.C[:, :, j]*d.Phi[:, :, d.tauIdx[i+1], d.tauIdx[j+1]]*d.A[:, :, d.tauIdx[i+1]-1]*d.xpts[:, d.tauIdx[i+1]])[1]
      # Htemp[j, i] = 0.5*(d.xpts[:, d.tauIdx[i+1]]'*d.A[:, :, d.tauIdx[i+1] - 1]'*d.Phi[:, :, d.tauIdx[i+1], d.tauIdx[j+1]]'*d.C[:, :, j]*d.xpts[:, d.tauIdx[j+1]]   +   d.xpts[:, d.tauIdx[j+1]]'*d.C[:, :, j]*d.Phi[:, :, d.tauIdx[i+1], d.tauIdx[j+1]]*d.A[:, :, d.tauIdx[i+1] - 1]*d.xpts[:, d.tauIdx[i+1]])[1]
    end
  end



# Try new indeces
  # for i = 1:d.N+1
  #   Htemp[i, i] = (d.xpts[:, d.tauIdx[i+1]]'*d.C[:, :, i]*d.A[:, :, d.tauIdx[i+1]-1]*d.xpts[:, d.tauIdx[i+1]])[1]
  # end
  #
  # for j = 2:d.N+1
  #   # for j = i+1:d.N+1
  #     for i = 1:j-1
  #     Htemp[j, i] = (d.xpts[:, d.tauIdx[j+1]]'*d.C[:, :, j]*d.Phi[:, :, d.tauIdx[i+1], d.tauIdx[j+1]]*d.A[:, :, d.tauIdx[i+1]-1]*d.xpts[:, d.tauIdx[i+1]])[1]
  #   end
  # end



  ### CONTINUE FROM HERE. Check Stuff. Change Simulations. Try Code.
  ### Maybe Change Linear Part as well and Try Code with that.

  H[:] = sigma*Htemp[d.IndTril]
end



# Constraints for Linear System

function MathProgBase.eval_g(d::linSTOev, g, x)

  # Constraint on sum of Switching Times
  g[1] = (d.gsum*x)[1]

  if d.ncons!=0  # Stagewise Constraints
    for k = 2:d.ngrid - 1  # iterate over grid
      # for l = 1:d.ncons    # Iterate over Stagewise constraints
        g[1+((k-2)*d.ncons + 1):1+(k-1)*d.ncons] = d.Ac*d.xpts[:, d.tgridIdx[k]]
      # end
    end
  end

  if d.nconsf!=0  # Final Stage Constraints
    g[end-d.nconsf+1:end] = d.Acf*d.xpts[:, end]
  end

 # print("\ng = $(g)\n\n")
  # # Old Constraints
  # if d.ncons!=0
  #   # Constraint on Last Stage
  #   g[2:end] = d.Ac*d.xpts[:, end]
  # end



  # Old one
  # g[:] = d.Ag*x

end


function MathProgBase.eval_jac_g(d::linSTOev, J, x)

Jac_temp = zeros(d.ncons*(d.ngrid - 2), d.N+1)

if d.ncons!=0  # There are stage constraints in the problem
    # Precompute matrices if necessary
    if d.prev_delta != x
        precompMatrices!(d, x) # Precompute Matrices and store them in d
        d.prev_delta[:] = x  # Update current tau
    end

    # println("Here!")
    # show(d.tauIdx)
    # show(d.tgridIdx)
    for k = 2:d.ngrid - 1   # Iterate over the grid constraints
      for l = 1:d.ncons    # Iterate over every element of the stagewise constraints
        for i = 1:d.N+1     # Iterate over the intervals
          # print("k = $(k), i = $(i), tauIdx[i+1]=$(d.tauIdx[i+1]), tgridIdx[k] = $(d.tgridIdx[k]) \n")
          if d.tauIdx[i+1] < d.tgridIdx[k]  # Change contraint only if not zero.
            Jac_temp[(k-2)*d.ncons+l, i] = (d.Ac[l, :]*d.Phi[:, :, d.tauIdx[i+1], d.tgridIdx[k]]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
          end
        end
      end
    end

end

Jac_tempf = zeros(d.nconsf, d.N+1)

if d.nconsf!=0  # There are  final stage constraints

  # Precompute matrices if necessary
  if d.prev_delta != x
      precompMatrices!(d, x) # Precompute Matrices and store them in d
      d.prev_delta[:] = x  # Update current tau
  end

  # Construct jacobian (Constraint on last stage)
  for l = 1:d.nconsf # Iterate over Constraints
    for i = 1:d.N+1 # Iterate over Variables
      Jac_tempf[l, i] = (d.Acf[l, :]*d.Phi[:, :, d.tauIdx[i+1], end]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
    end
  end

end

# Change jacobian if there are stage constraints
if (d.nconsf!=0) | (d.ncons!=0)
  # Compute Jacobian Matrix
  Jac =  [d.gsum; Jac_temp; Jac_tempf]
  d.Vg = Jac[:]
  # print("J = $(Jac)\n")
end

# show(d.nconsf); @printf("\n")
# show(d.ncons); @printf("\n")

  # Old Jacobian
  # if d.ncons!=0 # If there are any constraints in the problem
  #
  #   # Check if the matrices have already been precomputed
  #   if d.prev_delta != x
  #       precompMatrices!(d, x) # Precompute Matrices and store them in d
  #       d.prev_delta[:] = x  # Update current tau
  #   end
  #
  #
  #
  # # Compute Jacobian
  # Jac_temp = zeros(1+d.ncons, d.N+1)
  # Jac_temp[1,:] = d.gsum  # Constraint on the sum of Variables
  #
  # # Construct jacobian (Constraint on last stage)
  # for l = 1:d.ncons # Iterate over Constraints
  #   for i = 1:d.N+1 # Iterate over Variables
  #     Jac_temp[1+l, i] = (d.Ac[l, :]*d.Phi[:, :, d.tauIdx[i+1], end]*d.A[:, :, i]*d.xpts[:, d.tauIdx[i+1]])[1]
  #   end
  # end
  #
  # d.Vg = Jac_temp[:]
  #
  # end


  J[:] = d.Vg
end


function MathProgBase.jac_structure(d::linSTOev)
  # Check if the matrices have already been precomputed
  # if d.prev_delta != x
  #     precompMatrices!(d, x) # Precompute Matrices and store them in d
  #     d.prev_delta[:] = x  # Update current tau
  # end

  return d.Ig, d.Jg
end




# Constraints for nonlinear system (To Complete!)

function MathProgBase.eval_g(d::nlinSTOev, g, x)
  g[:] = d.Ag*x
end


function MathProgBase.eval_jac_g(d::nlinSTOev, J, x)
  J[:] = d.Vg
end


function MathProgBase.jac_structure(d::nlinSTOev)
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

  # display(ttemp)
  # @printf("tidxtemp = "); show(tidxtemp); @printf("\n")
  # display(tidxtemp)
  # display(tau)


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
function propagateDynamics!(d::linSTOev, tau::Array{Float64,1})
# function propagateDynamics!(d::linSTOev, tau::Array{Float64,1}, x::Array{Float64, 1})

  # Create merged and sorted time vector with grid and switching times
  d.tvec, d.tauIdx, d.tgridIdx = mergeSortFindIndex(d.tgrid, tau)

  # Get complete delta vector with all intervals
  d.deltacomplete = tau2delta(d.tvec[2:end-1], d.t0, d.tf)

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
function propagateDynamics!(d::nlinSTOev, tau::Array{Float64,1})

  # Fit switching times withing the grid
  d.tvec, d.tauIdx, d.tgridIdx = mergeSortFindIndex(d.tgrid, tau)

  # Get complete delta vector with all intervals
  d.deltacomplete = tau2delta(d.tvec[2:end-1], d.t0, d.tf)

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
