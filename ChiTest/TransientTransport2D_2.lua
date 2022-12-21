-- 1D Transient Transport test with Vacuum BC.
-- SDM: PWLD
-- Test:
num_procs = 2





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
N=4
L=100.0
xmin = 0.0
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

-- Define microscopic cross sections
xs_critical = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(xs_critical, CHI_XSFILE, "ChiTest/xs_inf_critical_1g.cxs")
xs_21cent = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(xs_21cent, CHI_XSFILE, "ChiTest/xs_inf_21cent_1g.cxs")

num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        EXISTING,xs_critical)

src={0.0}
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

function SwapXS(solver_handle, new_xs)
    chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
            EXISTING,new_xs)
    chiLBSInitializeMaterials(solver_handle)
end

--############################################### Setup Physics
phys1 = chiLBSCreateTransientSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
fac=1
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4*fac, 3*fac)
chiOptimizeAngularQuadratureForPolarSymmetry(pqaud, 4.0*math.pi)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,8)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,1000)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
--chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

--
----############################################### Set boundary conditions
--bsrc={}
--for g=1,num_groups do
--    bsrc[g] = 0.0
--end
--bsrc[1] = 1.0/2
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,LBSBoundaryTypes.REFLECTING);
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMAX,LBSBoundaryTypes.REFLECTING);
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMIN,LBSBoundaryTypes.REFLECTING);
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,YMAX,LBSBoundaryTypes.REFLECTING);
--
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)

chiLBKESSetProperty(phys1, "MAX_ITERATIONS", 100)
chiLBKESSetProperty(phys1, "TOLERANCE", 1.0e-8)

chiLBSSetProperty(phys1, USE_PRECURSORS, true)

--chiLBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
chiLBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
chiLBSSetProperty(phys1, VERBOSE_OUTER_ITERATIONS, true)


--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)

chiLBTSSetProperty(phys1, "TIMESTEP", 1e-3)
chiLBTSSetProperty(phys1, "VERBOSITY_LEVEL", 0)
chiLBTSSetProperty(phys1, "TIMESTEP_METHOD", "CRANK_NICHOLSON")

phys1name = chiSolverGetName(phys1);
initial_FR = chiLBSComputeFissionRate(phys1,"OLD")

--time = 0.0
--for k=1,2 do
--    --chiLBTSSetProperty(phys1, "INHIBIT_ADVANCE", true)
--    chiSolverStep(phys1)
--    FRf = chiLBSComputeFissionRate(phys1,"NEW")
--    FRi = chiLBSComputeFissionRate(phys1,"OLD")
--    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
--    time = chiLBTSGetProperty(phys1, "TIME")
--    period = dt/math.log(FRf/FRi)
--    chiLog(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
--            phys1name,k,time,dt,period,FRf/initial_FR))
--end

time = 0.0
time_stop = 1.0
k=0
swapped = false
while (time < time_stop) do
    k = k + 1
    chiSolverStep(phys1)
    FRf = chiLBSComputeFissionRate(phys1,"NEW")
    FRi = chiLBSComputeFissionRate(phys1,"OLD")
    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
    time = chiLBTSGetProperty(phys1, "TIME")
    period = dt/math.log(FRf/FRi)
    chiLog(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
            phys1name,k,time,dt,period,FRf/initial_FR))
    if (time >= 0.2 and not swapped) then
        SwapXS(phys1, xs_21cent)
        swapped = true
    end

    chiLBTSAdvanceTimeData(phys1)
end