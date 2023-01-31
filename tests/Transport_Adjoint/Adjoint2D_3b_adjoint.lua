-- 2D Transport test with point source Multigroup Adjoint generation
-- SDM: PWLD
-- Test: None
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
tmesh = chiMeshHandlerCreate()

nodes={}
N=60
L=5.0
ds=L/N
xmin=0.0
for i=0,N do
    nodes[i+1] = xmin + i*ds
end
mesh = chiMeshCreateUnpartitioned2DOrthoMesh(nodes,nodes)
chiVolumeMesherExecute();

----############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-1000,1000,0.0,0.8*L,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)



----############################################### Set Material IDs
vol0b = chiLogicalVolumeCreate(RPP,-0.166666+2.5,0.166666+2.5,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0b,0)

vol1b = chiLogicalVolumeCreate(RPP,-1+2.5,1+2.5,0.9*L,L,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1b,1)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 10
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,num_groups,0.01,0.01)
chiPhysicsMaterialSetProperty(materials[2],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,num_groups,0.1*20,0.8)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 1.0

--############################################### Setup Physics
solver_name = "LBAdjoint"
phys1 = chiAdjointSolverCreate(solver_name)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12, 2)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,500)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--############################################### Set boundary conditions

--############################################### Add point source
chiLBSAddPointSource(phys1, 1.25 - 0.5*ds, 1.5*ds, 0.0, src)

--############################################### Set solver properties
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)

--############################################### Create QOIs
tvol0 = chiLogicalVolumeCreate(RPP,2.3333,2.6666,4.16666,4.33333,-1000,1000)
tvol1 = chiLogicalVolumeCreate(RPP,0.5   ,0.8333,4.16666,4.33333,-1000,1000)

function ResponseFunction(cell_centroid, material_id)
    response = {}
    for g=0,9 do
        response[g+1] = 0.0
    end
    response[5+1] = 1.0

    return response
end

chiAdjointSolverAddResponseFunction(phys1,"QOI0",tvol0,"ResponseFunction")
chiAdjointSolverAddResponseFunction(phys1,"QOI1",tvol1,"ResponseFunction")
chiSolverSetBasicOption(phys1, "REFERENCE_RF", "QOI1")

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

chiLBSWriteFluxMoments(phys1, "Adjoint2D_3b_adjoint")



--############################################### Get field functions
ff_m0 = chiGetFieldFunctionHandleByName(solver_name.."_Flux_g0_m0")
ff_m1 = chiGetFieldFunctionHandleByName(solver_name.."_Flux_g0_m1")
ff_m2 = chiGetFieldFunctionHandleByName(solver_name.."_Flux_g0_m2")


--############################################### Slice plot

--############################################### Volume integrations












--############################################### Exports
if master_export == nil then
    chiExportMultiFieldFunctionToVTK({ff_m0, ff_m1, ff_m2},"ZPhi_"..solver_name)
    chiExportFieldFunctionToVTKG(ff_m0, "ZPhi_"..solver_name)
end
