-- 1D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.49903 and 7.18243e-4
num_procs = 3





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
N=100
L=30.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"ChiTest/xs_3_170.cxs")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"ChiTest/xs_3_170.cxs")

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
phys1 = chiLBSCreateSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE,40)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,62)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,8)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
--chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

gs1 = chiLBSCreateGroupset(phys1)
cur_gs = gs1
chiLBSGroupsetAddGroups(phys1,cur_gs,63,167)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,8)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
--chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

--############################################### Set boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/2
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,5)

--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--############################################### Line plot
--Testing consolidated interpolation
cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,0.0,0.0001+xmin)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,0.0, 29.999+xmin)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 50)

for k=165,165 do
    chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[k])
end

chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)

--############################################### Volume integrations
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value1=%.5f", maxval))

ffi2 = chiFFInterpolationCreate(VOLUME)
curffi = ffi2
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[160])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

--############################################### Exports
if (master_export == nil) then
    chiFFInterpolationExportPython(cline)
end

--############################################### Plots
if (chi_location_id == 0 and master_export == nil) then
    local handle = io.popen("python3 ZLFFI00.py")
end
