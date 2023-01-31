#include "ChiLua/chi_lua.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "chi_runtime.h"
#include "chi_log.h"
;

#include "chi_runtime.h"
#include "ChiTimer/chi_timer.h"


#include "ChiConsole/chi_console.h"

#include <iomanip>
#include <iostream>
#include <unistd.h>

//#############################################################################
/** Executes the volume meshing pipeline.

\ingroup LuaVolumeMesher
\author Jan*/
int chiVolumeMesherExecute(lua_State *L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //Get memory before
  chi_objects::CSTMemory mem_before = chi::console.GetMemoryUsage();

  if (cur_hndlr.volume_mesher == nullptr)
  {
    chi::log.LogAllError()
      << __FUNCTION__ << ": called without a volume mesher set. Make a "
                         "call to chiVolumeMesherCreate.";
    chi::Exit(EXIT_FAILURE);
  }

  cur_hndlr.volume_mesher->Execute();

  //Get memory usage
  chi_objects::CSTMemory mem_after = chi::console.GetMemoryUsage();

  std::stringstream mem_string;
  mem_string
  << " Memory used = " << std::setprecision(3)
  << mem_after.memory_mbytes - mem_before.memory_mbytes
  << " MB\n"
  << "Total process memory used after meshing "
  << mem_after.memory_mbytes
  << " MB";

  chi::log.Log()
    << chi::program_timer.GetTimeString()
    << " chiVolumeMesherExecute: Volume meshing completed."
    << mem_string.str()
    << std::endl;

  return 0;
}