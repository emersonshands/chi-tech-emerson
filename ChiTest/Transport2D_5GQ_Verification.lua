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
sn = 4
method = 1

scatterOrder = 3 --2*(sn-1)
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2,2)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad,4.0*math.pi)
pq2 = chiCreateProductQuadratureOperator(pquad,method,sn)

tab = chiGetProductQuadrature(pq2)
chiLog(LOG_0, "Checking Values of Quadrature")
weights = {}
for pl=1,rawlen(tab) do
    chiLog(LOG_0, "Direction " .. tostring(pl))
    chiLog(LOG_0, "Weight " .. tostring(tab[pl].weight))
    chiLog(LOG_0, "Polar " .. tostring(tab[pl].polar))
    chiLog(LOG_0, "XCos " .. tostring(math.sin(tab[pl].polar)*math.cos(tab[pl].azimuthal)*tab[pl].weight))
    weights[#weights + 1] = tab[pl].weight
    chiLog(LOG_0, "Azimu " .. tostring(tab[pl].azimuthal))
end

--scatterOrder = 2
--pquad = chiCreateAngularQuadratureTriangle(method,sn)
--
--tab = chiGetTriangleQuadrature(pquad)
--chiLog(LOG_0, "Checking Values of Quadrature")
--weights = {}
--for pl=1,rawlen(tab) do
--    chiLog(LOG_0, "Direction " .. tostring(pl))
--    chiLog(LOG_0, "Weight " .. tostring(tab[pl].weight))
--    weights[#weights + 1] = tab[pl].weight
--    chiLog(LOG_0, "Polar " .. tostring(tab[pl].polar))
--    chiLog(LOG_0, "Azimu " .. tostring(tab[pl].azimuthal))
--end

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pq2)
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
    sum = 0
    for ni=1,num_angles do
        --print("Sum " .. tostring(sum))
        omega = quadrature_angle_vectors[ni]
        phi_theta = quadrature_phi_theta_angles[ni]
        indexQ = quadrature_angle_indices[ni]
        --print("Current Index " .. tostring(indexQ))
        weightCurrent = weights[indexQ+1]
        --print("Current weight " .. tostring(weightCurrent) .. " Current index " .. tostring(indexQ))
        --print("Current Phi " .. tostring(phi_theta.phi*180/math.pi))
        for gi=1,num_groups do
            g = group_indices[gi]

            value = 0.0
            if ( omega.x > 0 and omega.y > 0 and omega.z >0) then
                value = 1.0
                sum = sum + omega.x*weightCurrent
            end
            --print("Used Value " .. tostring(value) .. " Dot product " .. tostring(dot))
            dof_count = dof_count + 1
            psi[dof_count] = value
            --print("Sum " .. tostring(sum))
        end
    end
    for angl=1,dof_count do
        psi[angl] = psi[angl]/sum/10
        --print("psi value " .. tostring(psi[angl]))
    end
    return psi
end

----############################EXTRA BC##################################
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
----############################EXTRA BC##################################
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
chiLBSSetProperty(phys1,SCATTERING_ORDER,scatterOrder)
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
if not master_export == nil then
    chiFFInterpolationExportPython(slice2)
end

leakage0 = chiLBSComputeLeakage(phys1, gs0, 0)
leakage1 = chiLBSComputeLeakage(phys1, gs0, 1)
leakage2 = chiLBSComputeLeakage(phys1, gs0, 2)
leakage3 = chiLBSComputeLeakage(phys1, gs0, 3)
balance = chiLBSComputeBalance(phys1)

--chiPrintD2M(pquad)
--chiPrintM2D(pquad)
--chiCheckIdentity(pquad)

chiLog(LOG_0,"XMax")
chiLog(LOG_0,tostring(leakage0[1]))
chiLog(LOG_0,"XMin")
chiLog(LOG_0,tostring(leakage1[1]))
chiLog(LOG_0,"YMax")
chiLog(LOG_0,tostring(leakage2[1]))
chiLog(LOG_0,"Ymin")
chiLog(LOG_0, tostring(leakage3[1]))

--############################################### Plots
if (chi_location_id == 0 and not master_export == nil) then
    local handle = io.popen("python3 ZPFFI00.py")
end

