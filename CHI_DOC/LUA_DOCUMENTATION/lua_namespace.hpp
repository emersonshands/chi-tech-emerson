namespace CHI_LUA 
 {
int chiPieExportPin(int pinNumber);
int chiPieSetPinMode(int pinNumber, int mode);
int chiPieSetPinValue(int pinNumber, int mode);
int chiPieGetPinValue(int pinNumber);
int chiPieInitSPI();
int chiPieReadSPIChannel(int channelNumber);
int chiPieSetSPIBuffer(int channelNumber, bool bufferFlag);
int chiPieGetSPIBuffer(int channelNumber, int bufferPos);
int chiPieInitializeSerial(int baudrate);
int chiPieSerialWrite(char message);
int chiPieSerialRead();
int chiThermoSetComponentProperty(Handle sysHndle, Handle compHndle, Property propCode);
int chiThermoCreateVolumeFromCoordinates(int systemHandle, Table point1, Table point2);
int chiThermoCreateSystem();
int chiThermoConnectTwoComponents(Handle systemHandle, Component leftComponent, Single sjunc, Component rigtComponent, 0=end-begin, mode);
int chiThermoInitialize(int systemHandle);
int chiThermoCreateBC(int systemHandle);
int chiThermoGetComponentProperty(Handle sysHndle, Handle compHndle, Property propCode);
int chiThermoCreateSJunction(Handle systemHandle);
int chiCreateProductQuadrature(int QuadratureType, int Np, int Na);
int chiCreateQuadrature(int QuadratureType, int NumberOfPoints);
int chiLegendre(int N, double x);
int chiLegendreDerivative(int N, double x);
int chiYlm(ell \param, m \param, theta \param, varphi \param);
int chiEdgeLoopSplitByAngle(int LoopCollectionHandle, int LoopHandle, double Angle);
int chiSurfaceMeshCreate();
int chiSurfaceMeshImportFromOBJFile(int SurfaceHandle, char* FileName, bool polyflag);
int chiSurfaceMeshExtractOpenEdgesToObj(int SurfaceHandle, char FileName);
int chiSurfaceMeshCheckCycles(int SurfaceHandle, int NumAngles);
int chiSurfaceMeshSplitByPatch(int SurfaceHandle);
int chiSurfaceMeshGetEdgeLoops(int SurfaceHandle);
int chiSurfaceMeshGetEdgeLoopsPoly(int SurfaceHandle);
int chiLogicalVolumeCreate(int TypeIndex, varying Values);
int chiRegionCreate();
int chiRegionGetBoundarySurfaceMesh(int RegionHandle, int BoundaryNumber, int ContinuumNumber);
int chiRegionExportMeshToPython(int RegionHandle, char FileName, bool ExportTemplate);
int chiRegionExportMeshToObj(int RegionHandle, char FileName, bool ExportByMaterial);
int chiRegionExportMeshToVTK(int RegionHandle, char FileName);
int chiRegionAddSurfaceBoundary(int RegionHandle, int SurfaceHandle);
int chiRegionAddLineBoundary(int RegionHandle, int LineMeshHandle);
int chiMeshHandlerCreate();
int chiMeshHandlerSetCurrent(int HandlerHandler);
int chiMeshHandlerGetSurfaceFromCollection(int CollectionHandle, int SurfaceIndex);
int chiLineMeshCreateFromArray(LuaTable Table);
int chiLineMeshCreateFromLoop(int LoopCollectionHandle, int LoopHandle);
int chiSurfaceMesherExportToObj();
int chiSurfaceMesherSetProperty(int PropertyNumber, varying PropertyValue);
int chiSurfaceMesherExecute(bool ExportLoadBalance);
int chiSurfaceMesherCreate(int Type);
int chiVolumeMesherCreate(int Type);
int chiVolumeMesherExecute();
int chiVolumeMesherSetProperty(int PropertyIndex, varying PropertyValue);
int chiFFInterpolationInitialize(int FFIHandle);
int chiFFInterpolationExecute(int FFIHandle);
int chiFFInterpolationGetValue(int FFIHandle);
int chiFFInterpolationCreate(int FFITypeIndex);
int chiFFInterpolationSetProperty(int FFIHandle, int PropertyIndex);
int chiFFInterpolationExportPython(int FFIHandle, char BaseName);
int chiMPIReceiveCellsets();
int chiMPIBroadcastCellsets();
int chiMPIBarrier();
int chiLogSetVerbosity(int int_level);
int chiLog(int LogType);
int chiPhysicsAddMaterial(char Name);
int chiExportFieldFunctionToVTK(int FFHandle, char BaseName);
int chiExportFieldFunctionToVTKG(int FFHandle, char BaseName);
int chiSolverExecute();
int chiSolverAddFieldFunction(int SolverHandle, char Name);
int chiPhysicsMaterialAddProperty(int MaterialHandle, int PropertyIndex);
int chiPhysicsMaterialSetProperty(int MaterialHandle, int PropertyIndex, int OperationIndex, varying Information);
int chiGetFieldFunctionList(SolverHandle \param);
int chiSolverAddRegion(int SolverHandle, int RegionHandle);
int chiMonteCarlonCreateSolver();
int chiMonteCarlonCreateSource(int SolverHandle, int SourceType);
int chiMonteCarlonInitialize(int SolverHandle);
int chiMonteCarlonSetProperty(int SolverHandle, int PropertyIndex);
int chiMonteCarlonExecute(int SolverHandle);
int chiDiffusionSetProperty(int SolverHandle, int PropertyIndex, varying Values);
int chiDiffusionCreateSolver();
int chiDiffusionInitialize(int SolverHandle);
int chiDiffusionExecute(int SolverHandle, int SolverHandle);
int chiNPTInitialize(int SolverIndex);
int chiNPTransportCreateSolver();
int chiNPTGetFieldFunctionList(int SolverIndex);
int chiNPTGetScalarFieldFunctionList(int SolverIndex);
int chiNPTCreateMaterial(int SolverIndex);
int chiNPTMaterialSetPureAbsorber(int SolverIndex, int MaterialIndex, int NumberOfGroups, int ScatOrder, double Sigma_t=1.0);
int chiNPTSetProperty(int SolverIndex, int PropertyIndex);
int chiNPTCreateGroupset(int SolverIndex);
int chiNPTCreateGroup(int SolverIndex);
int chiNPTGroupsetAddGroups(int SolverIndex, int GroupsetIndex, int FromIndex, int ToIndex);
int chiNPTGroupsetSetQuadrature(int SolverIndex, int GroupsetIndex, int QuadratureIndex);
int chiNPTGroupsetSetAngleAggDiv(int SolverIndex, int GroupsetIndex, int NumDiv);
int chiNPTGroupsetSetGroupSubsets(int SolverIndex, int GroupsetIndex, int NumDiv);
int chiNPTGroupsetSetIterativeMethod(int SolverIndex, int GroupsetIndex, int IterativeMethod);
int chiNPTGroupsetSetResidualTolerance(int SolverIndex, int GroupsetIndex, float ResidualTol);
int chiNPTGroupsetSetMaxIterations(int SolverIndex, int GroupsetIndex, int Numiter);
int chiNPTGroupsetSetGMRESRestartIntvl(int SolverIndex, int GroupsetIndex, int Intvl);
int chiNPTGroupsetSetWGDSA(int SolverIndex, int GroupsetIndex, int MaxIters, float ResTol, bool Verbose, char PETSCString);
int chiNPTGroupsetSetTGDSA(int SolverIndex, int GroupsetIndex, int MaxIters, float ResTol, bool Verbose, char PETSCString);
int chiNPTExecute(int SolverIndex);
}