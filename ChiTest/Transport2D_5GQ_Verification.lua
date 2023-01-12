-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
num_procs = 1


--############################################### Check num_procs
--if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
--    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
--            "Expected "..tostring(num_procs)..
--            ". Pass check_num_procs=false to override if possible.")
--    os.exit(false)
--end
--############################################### Setup mesh
--chiMeshHandlerCreate()
--chiUnpartitionedMeshFromWavefrontOBJ("ChiResources/TestObjects/TriangleMesh2x2.obj")
--region1 = chiRegionCreate()
--chiRegionAddEmptyBoundary(region1)
--chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
--chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED)
--chiSurfaceMesherExecute()
--chiVolumeMesherExecute()
chiMeshHandlerCreate()
mesh={}
N=65
L=10.0
xmin = 0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
--chiVolumeMesherSetProperty(PARTITION_TYPE,KBA_STYLE_XYZ)
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
        CHI_XSFILE,"ChiTest/xs_quad_test_GQ.cxs")

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
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2,2)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad,4.0*math.pi)
--pquad = chiCreateProductQuadratureOperator(pbase,3,4)
--pquad = chiCreateAngularQuadratureTriangle(1,4,1)
--chiPrintD2M(pquad)
--chiPrintM2D(pquad)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,2.606e-10)
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


function luaBoundaryFunctionLeft(cell_global_id,
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
    anglePass = 0

    normalVal = 1.0
    for ni=1,num_angles do
        omega = quadrature_angle_vectors[ni]
        phi_theta = quadrature_phi_theta_angles[ni]
        for gi=1,num_groups do
            g = group_indices[gi]

            value = normalVal
            --if (phi_theta.phi<1.57079632679 or phi_theta.phi>4.71238898038) then
            --if (phi_theta.phi<1.57079632679 and phi_theta.phi>0) then
            --    value = normalVal
            --end

            dof_count = dof_count + 1
            psi[dof_count] = value
        end
    end
    return psi
end
function luaBoundaryFunctionRight(cell_global_id,
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
    anglePass = 0

    normalVal = 1.0/4.0/3.14159265359/1.0
    for ni=1,num_angles do
        omega = quadrature_angle_vectors[ni]
        phi_theta = quadrature_phi_theta_angles[ni]
        for gi=1,num_groups do
            g = group_indices[gi]

            value = 0.0
            if (phi_theta.phi>1.57079 and phi_theta.phi<4.71238) then
                value = normalVal
            end

            dof_count = dof_count + 1
            psi[dof_count] = value
        end
    end
    return psi
end
function luaBoundaryFunctionBottom(cell_global_id,
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
    anglePass = 0

    normalVal = 1.0/4.0/3.14159265359/1.0
    for ni=1,num_angles do
        omega = quadrature_angle_vectors[ni]
        phi_theta = quadrature_phi_theta_angles[ni]
        for gi=1,num_groups do
            g = group_indices[gi]

            value = 0.0
            if (phi_theta.phi>0 and phi_theta.phi<3.14159265359) then
                value = normalVal
            end

            dof_count = dof_count + 1
            psi[dof_count] = value
        end
    end
    return psi
end
function luaBoundaryFunctionTop(cell_global_id,
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
    anglePass = 0

    normalVal = 1.0/4.0/3.14159265359/1.0
    for ni=1,num_angles do
        omega = quadrature_angle_vectors[ni]
        phi_theta = quadrature_phi_theta_angles[ni]
        for gi=1,num_groups do
            g = group_indices[gi]

            value = 0.0
            if (phi_theta.phi>3.14159265359 and phi_theta.phi<6.28318530718) then
                value = normalVal
            end

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
                        "luaBoundaryFunctionLeft");
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMAX,
--        LBSBoundaryTypes.INCIDENT_ANISTROPIC_HETEROGENOUS,
--        "luaBoundaryFunctionRight");
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMAX,
--        LBSBoundaryTypes.INCIDENT_ANISTROPIC_HETEROGENOUS,
--        "luaBoundaryFunctionTop");
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMIN,
--        LBSBoundaryTypes.INCIDENT_ANISTROPIC_HETEROGENOUS,
--        "luaBoundaryFunctionBottom");

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)
chiLBSSetProperty(phys1,SAVE_ANGULAR_FLUX,true)

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

--############################################### Exports
if master_export == nil then
    chiFFInterpolationExportPython(slice2)
end

leakage0 = chiLBSComputeLeakage(phys1, gs0, 0)
leakage1 = chiLBSComputeLeakage(phys1, gs0, 1)
leakage2 = chiLBSComputeLeakage(phys1, gs0, 2)
leakage3 = chiLBSComputeLeakage(phys1, gs0, 3)
balance = chiLBSComputeBalance(phys1)

print("XMax")
chiLog(LOG_0,tostring(leakage0[1]))
print("XMin")
chiLog(LOG_0,tostring(leakage1[1]))
print("YMax")
chiLog(LOG_0,tostring(leakage2[1]))
print("YMin")
chiLog(LOG_0, tostring(leakage3[1]))

--############################################### Plots
if (chi_location_id == 0 and master_export == nil) then
    local handle = io.popen("python3 ZPFFI00.py")
end





