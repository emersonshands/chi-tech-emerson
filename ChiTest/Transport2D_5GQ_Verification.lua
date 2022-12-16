-- 2D Transport to verify GQ Paper on M2D and D2M methods.
num_procs = 10

--############################################### Check num_procs
--if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
--    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
--            "Expected "..tostring(num_procs)..
--            ". Pass check_num_procs=false to override if possible.")
--    os.exit(false)
--end

--############################################### Setup mesh
chiMeshHandlerCreate()
nodes={}
N=65
lent = 10
ds=lent/N
for i=0,N do
    nodes[i+1] = i*ds
end
surf_mesh,region1 = chiMeshCreateUnpartitioned2DOrthoMesh(nodes,nodes)

chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000E6,1000E6,-1000E6,1000E6,-1000E6,1000E6)
--Get the right side volume
vol1 = chiLogicalVolumeCreate(RPP,lent-ds,lent,0,lent,-1000E6,1000E6)
--Get the top side volume
vol2 = chiLogicalVolumeCreate(RPP,0,lent,lent-ds,lent,-1000E6,1000E6)
--Get the bottom volume
vol3 = chiLogicalVolumeCreate(RPP,0,lent,0,ds,-1000E6,1000E6)
--Get the left volume
vol4 = chiLogicalVolumeCreate(RPP,-ds,0,0,lent,-1000E6,1000E6)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Scatter Material");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"ChiTest/xs_quad_test_GQ.cxs")

--############################################### Setup Physics

phys1 = chiLBSCreateSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquadOp = chiCreateAngularQuadratureTriangle(4,1)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquadOp)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--############################################### Set boundary conditions

--This needs to be a "unit incident half-range current"
bsrc={}
for g=1,num_groups do
    bsrc[g] = 1.0;
end
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);


chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)

--############################################### Initialize and Execute Solver
chiLBSInitialize(phys1)
chiLBSExecute(phys1)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--############################################### Volume integrations
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1

--chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])
--
--chiFFInterpolationSetProperty(curffi,LINE_FIRSTPOINT,lent,lent,0)
--chiFFInterpolationSetProperty(curffi,LINE_SECONDPOINT,lent,0,0)
--chiFFInterpolationSetProperty(curffi,LINE_NUMBEROFPOINTS,65)
--chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol1)
--

chiFFInterpolationSetProperty(curffi,OPERATION,OP_SUM)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol1)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])


chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
RightSideSum = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("RightSum=%.5e", RightSideSum))
--Desired value for this right side sum is [0.0199, 0.0182]
--This gets 0.0801 currently doing the volume integration
--A way to get the boundary leakage needs to be found.

--############################################### Exports

chiExportFieldFunctionToVTKG(fflist[1],"Phi2D_GQ","Phi")


