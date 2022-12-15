-- 3D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=5.27450e-01 and 3.76339e-04
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
N=200
lent = 1E5
ds=lent/N
for i=0,N do
    nodes[i+1] = i*ds
end
surf_mesh,region1 = chiMeshCreateUnpartitioned2DOrthoMesh(nodes,nodes)

chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000E6,1000E6,-1000E6,1000E6,-1000E6,1000E6)
vol1 = chiLogicalVolumeCreate(RPP,lent/4.0,lent*3.0/4.0,lent/4.0,lent*3.0/4.0,-1000E6,1000E6)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)



num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"ChiTest/xs_quad_test_scatter.cxs")

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 1.0

chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics

phys1 = chiLBSCreateSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
--chiLog(LOG_0,"Creating GLC quadratures")
----pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,1,1)
--chiLog(LOG_0,"Altering Quadrature")
-- pquadOp2 = chiCreateProductQuadratureOperator(pquad,3,2,1)
pquadOp = chiCreateAngularQuadratureTriangle(8,3)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
--chiLBSAddPointSource(phys1,lent/2.0,lent/2.0,0.0,{1.0/4.0/math.pi})
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquadOp)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
if (master_export == nil) then
    --chiLBSGroupsetSetEnableSweepLog(phys1,cur_gs,true)
end
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--############################################### Set boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 1.0/4.0/math.pi;
end
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMAX,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMAX,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMAX,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);


chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)

--############################################### Initialize and Execute Solver
chiLBSInitialize(phys1)
chiLBSExecute(phys1)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)
--############################################### Slice plot
slices = {}
for k=1,count do
    slices[k] = chiFFInterpolationCreate(SLICE)
    --chiFFInterpolationSetProperty(slices[k],SLICE_POINT,0.0,0.0,0.8001)
    chiFFInterpolationSetProperty(slices[k],ADD_FIELDFUNCTION,fflist[k])
    --chiFFInterpolationSetProperty(slices[k],SLICE_TANGENT,0.393,1.0-0.393,0)
    --chiFFInterpolationSetProperty(slices[k],SLICE_NORMAL,-(1.0-0.393),-0.393,0.0)
    chiFFInterpolationSetProperty(slices[k],SLICE_BINORM,0.0,0.0,1.0)
    chiFFInterpolationInitialize(slices[k])
    chiFFInterpolationExecute(slices[k])
    chiFFInterpolationExportPython(slices[k])
end

--############################################### Volume integrations
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol1)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value1=%.5e", maxval))

ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[20])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

--############################################### Quadrature print
--chiPrintM2D(pquadOp);
--chiPrintD2M(pquadOp);
--############################################### Exports
chiExportFieldFunctionToVTKG(fflist[1],"Phi2D_Scatter","Phi")

--############################################### Plots
if (chi_location_id == 0 and master_export == nil) then

    --os.execute("python ZPFFI00.py")
    ----os.execute("python ZPFFI11.py")
    --local handle = io.popen("python ZPFFI00.py")
    print("Execution completed")
end

