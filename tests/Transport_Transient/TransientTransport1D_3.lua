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
N=2000
L=100.0
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)

vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-L/4,L/4)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,1)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Strong fuel");
materials[2] = chiPhysicsAddMaterial("Weak fuel");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)

-- Define microscopic cross sections
xs_strong_fuel_micro = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(xs_strong_fuel_micro, CHI_XSFILE, "tests/Transport_Transient/xs_inf_k1_6_1g.cxs")
xs_weak_fuelA_micro = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(xs_weak_fuelA_micro, CHI_XSFILE, "tests/Transport_Transient/xs_inf_critical_1g.cxs")
xs_weak_fuelB_micro = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(xs_weak_fuelB_micro, CHI_XSFILE, "tests/Transport_Transient/xs_inf_weak_1g.cxs")

atom_density = 0.056559
xs_strong_fuel = chiPhysicsTransportXSMakeCombined({{xs_strong_fuel_micro, atom_density}}) --critical
xs_weak_fuelA  = chiPhysicsTransportXSMakeCombined({{ xs_weak_fuelA_micro, atom_density}}) --critical
xs_weak_fuelB  = chiPhysicsTransportXSMakeCombined({{xs_weak_fuelB_micro, atom_density}})   --critical

num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        EXISTING,xs_strong_fuel)
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        EXISTING,xs_weak_fuelA)

src={0.0}
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

function SwapXS(solver_handle, new_xs)
    chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
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
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE,16)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,8)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES)
--chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_RICHARDSON)
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
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.REFLECTING);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMAX,LBSBoundaryTypes.REFLECTING);
--
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)

chiLBKESSetProperty(phys1, "MAX_ITERATIONS", 1000)
chiLBKESSetProperty(phys1, "TOLERANCE", 1.0e-8)

chiLBSSetProperty(phys1, USE_PRECURSORS, true)

--chiLBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
chiLBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
chiLBSSetProperty(phys1, VERBOSE_OUTER_ITERATIONS, true)


--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)

chiLBTSSetProperty(phys1, "TIMESTEP", 1e-3*100)
chiLBTSSetProperty(phys1, "TIMESTOP", 1.0*100)
chiLBTSSetProperty(phys1, "MAX_TIMESTEPS", -1)
chiLBTSSetProperty(phys1, "VERBOSITY_LEVEL", 0)
chiLBTSSetProperty(phys1, "TIMESTEP_METHOD", "CRANK_NICHOLSON")

phys1name = chiSolverGetName(phys1);
initial_FR = chiLBSComputeFissionRate(phys1,"OLD")


--time = 0.0
--psi_t = psi_0
--for k=1,10 do
--    chiSolverStep(phys1)
--
--    FRf = chiLBSComputeFissionRate(phys1,"NEW") --time+dt
--    FRi = chiLBSComputeFissionRate(phys1,"OLD") --time
--    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
--    time = chiLBTSGetProperty(phys1, "TIME")
--    new_time = time+dt
--
--    period = dt/math.log(FRf/FRi)
--    chiLog(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
--            phys1name,k,time,dt,period,FRf/initial_FR))
--
--    if (not timestep_rejected) then
--        chiLBTSAdvanceTimeValues()
--    end
--    chiLBTSAdvanceTimeData(phys1)
--end

--time = 0.0
--time_stop = 1.0*100
--k=0
--swapped = false
--timestep_rejected = false
--
--tolA = 10.0
--while (time < time_stop) do
--    chiSolverStep(phys1)
--    FRf = chiLBSComputeFissionRate(phys1,"NEW")
--    FRi = chiLBSComputeFissionRate(phys1,"OLD")
--    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
--    time = chiLBTSGetProperty(phys1, "TIME")
--    period = dt/math.log(FRf/FRi)
--
--    if (time >= 0.2 and not swapped) then
--        SwapXS(phys1, xs_weak_fuelB)
--        swapped = true
--    end
--
--    if (not timestep_rejected) then
--        chiLBTSAdvanceTimeData(phys1)
--        k = k + 1
--        chiLog(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
--                phys1name,k,time,dt,period,FRf/initial_FR))
--    else
--        timestep_rejected = false
--    end
--end
swapped=false
function MyCallBack()
    FRf = chiLBSComputeFissionRate(phys1,"NEW")
    FRi = chiLBSComputeFissionRate(phys1,"OLD")
    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
    time = chiLBTSGetProperty(phys1, "TIME")
    period = dt/math.log(FRf/FRi)

    if (time >= 0.2 and not swapped) then
        SwapXS(phys1, xs_weak_fuelB)
        swapped = true
    end
    chiLog(LOG_0,string.format("%s time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
                phys1name,time,dt,period,FRf/initial_FR))
end
chiLBTSSetProperty(phys1, "CALLBACK", "MyCallBack")
chiSolverExecute(phys1)