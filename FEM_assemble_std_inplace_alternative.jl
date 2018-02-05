# # ---------------------------------------------------------------------------------------------
# # FEM assembling
# # ---------------------------------------------------------------------------------------------

# ########################################
# ########################################
# # --------------------------
# # -------   Mass   -------
# # --------------------------

# function assemble_mass!(M  :: SparseMatrixCSC{Float64,Int64},
#                         m :: Mesh.TriangleMesh.TriMesh,
#                         d :: FEM.AbstractDof,
#                         r :: FEM.RefEl_Pk,
#                         q :: Quad.AbstractQuad,
#                         par :: Parameter.Parameter_FEM,
#                         p :: Problem.AbstractProblem)

#     n = r.n_node
    
#     # fixed quantities for mesh
#     weight_elem = FEM.map_ref_point_grad_det(m, q.point, 1:d.n_elem)

#     Phi = FEM.eval(r, q.point)
#     Phi_test = FEM.eval(r, q.point) * diagm(q.weight)
    
#     m_loc = Array{Float64,2}(n,n)
#     for k = 1:m.n_cell
#         for j = 1:length(q.weight)          
#             for ii in 1:n
#                 for ii_test in 1:n
#                     m_loc[ii_test,ii] = Phi_test[ii_test,j] * weight_elem[j,k] * Phi[ii,j]
#                 end
#             end

#             for l in 1:9
#                 M[d.ind_test[9*(k-1)+l],d.ind[9*(k-1)+l]] += m_loc[l]
#             end
#         end
#     end

# return nothing
# end
# ########################################
# ########################################



# ########################################
# ########################################
# # ---------------------------------
# # -------   Advection   -------
# # ---------------------------------

# ### This is buggy
# function assemble_advection!(A  :: SparseMatrixCSC{Float64,Int64},
#                m :: Mesh.TriangleMesh.TriMesh,
#                d :: FEM.AbstractDof,
#                r :: FEM.RefEl_Pk,
#                q :: Quad.AbstractQuad,
#                par :: Parameter.Parameter_FEM,
#                p :: Problem.AbstractProblem)

#     n = r.n_node
    
#     # fixed quantities for mesh
#     weight_elem = FEM.map_ref_point_grad_det(m, q.point, 1:d.n_elem)
#     DF = FEM.map_ref_point_grad_inv(m, q.point, 1:d.n_elem);

#     x = FEM.map_ref_point(m, q.point, 1:d.n_elem)
#     velocity = Problem.velocity(p, 1.0, x)

#     DPhi = FEM.eval_grad(r, q.point)
#     Phi_test = FEM.eval(r, q.point) * diagm(q.weight)

#     error("This function has a bug")
    
#     a_loc = Array{Float64,2}(n,n)
#     for k = 1:m.n_cell
#         for j = 1:length(q.weight)
#             for ii_test in 1:n
#                 for ii in 1:n
#                     a_loc[ii_test,ii] = Phi_test[ii_test,j] * weight_elem[j,k] * (DPhi[ii,j,1]*velocity[j,1,k] + DPhi[ii,j,2]*velocity[j,2,k])
#                 end
#             end   
#             #a_loc[:,:] = Phi_test[:,j] * weight_elem[j,k] * (DPhi[:,j,:]*velocity[j,:,k])'
#             for l in 1:9
#                 A[d.ind_test[9*(k-1)+l],d.ind[9*(k-1)+l]] += a_loc[l]
#             end
#         end
#     end
    
#     return nothing
# end
# ########################################
# ########################################



# ########################################
# ########################################
# # -------------------------------
# # -------   Diffusion   -------
# # -------------------------------


# ########################################
# ########################################
