dofile("mesh_3dhex.lua")

chi_unit_tests.chi_math_SDM_Test01_Continuous
({
  sdm_type = "LagrangeC",
  --export_vtk = true
});