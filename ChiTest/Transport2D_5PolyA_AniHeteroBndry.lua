-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.50758 and 2.52527e-04
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=40
L=10.0
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)
--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"ChiTest/xs_air50RH.cxs")

src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
phys1 = chiLBSCreateSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12, 2)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaud, 4.0*math.pi)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,0)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,2)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--gs1 = chiLBSCreateGroupset(phys1)
--cur_gs = gs1
--chiLBSGroupsetAddGroups(phys1,cur_gs,63,167)
--chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
--chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
--chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,2)
--chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
--chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
--chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
--chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
----chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
----chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

--############################################### Set boundary conditions
--int cell_global_id
--int material_id

--VecXYZ location (.x .y and .z)
--VecXYZ normal

--array<int>      quadrature_angle_indices
--array<VecXYZ>   quadrature_angle_vectors
--array<PhiTheta> quadrature_phi_theta_angles (PhiTheta.phi and PhiTheta.theta)
--array<int>      group_indices

function luaBoundaryFunctionA(cell_global_id,
                              material_id,
                              location,
                              normal,
                              quadrature_angle_indices,
                              quadrature_angle_vectors,
                              quadrature_phi_theta_angles,
                              group_indices)
    num_angles = rawlen(quadrature_angle_vectors)
    num_groups = rawlen(group_indices)
    psi = {}
    dof_count = 0

    for ni=1,num_angles do
        omega = quadrature_angle_vectors[ni]
        phi_theta = quadrature_phi_theta_angles[ni]
        for gi=1,num_groups do
            g = group_indices[gi]

            value = 1.0

            dof_count = dof_count + 1
            psi[dof_count] = value
        end
    end

    return psi
end
--############################################### Set boundary conditions


bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,
--        LBSBoundaryTypes.INCIDENT_ISOTROPIC, bsrc);
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,
                        LBSBoundaryTypes.INCIDENT_ANISTROPIC_HETEROGENOUS,
                        "luaBoundaryFunctionA");

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--############################################### Slice plot
slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)

----############################################### Volume integrations
--ffi1 = chiFFInterpolationCreate(VOLUME)
--curffi = ffi1
--chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
--chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
--chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])
--
--chiFFInterpolationInitialize(curffi)
--chiFFInterpolationExecute(curffi)
--maxval = chiFFInterpolationGetValue(curffi)
--
--chiLog(LOG_0,string.format("Max-value1=%.5f", maxval))

----############################################### Volume integrations
--ffi1 = chiFFInterpolationCreate(VOLUME)
--curffi = ffi1
--chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
--chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
--chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[160])
--
--chiFFInterpolationInitialize(curffi)
--chiFFInterpolationExecute(curffi)
--maxval = chiFFInterpolationGetValue(curffi)
--
--chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

--############################################### Exports
if master_export == nil then
    chiFFInterpolationExportPython(slice2)
end

--############################################### Plots
if (chi_location_id == 0 and master_export == nil) then
    local handle = io.popen("python ZPFFI00.py")
end