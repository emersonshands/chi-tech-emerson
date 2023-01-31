#ifndef CHI_SWEEPSCHEDULER_H
#define CHI_SWEEPSCHEDULER_H

#include "ChiMesh/SweepUtilities/AngleAggregation/angleaggregation.h"
#include "ChiMesh/SweepUtilities/sweepchunk_base.h"


namespace chi_mesh::sweep_management
{
  enum class SchedulingAlgorithm
  {
    FIRST_IN_FIRST_OUT = 1, ///< FIFO
    DEPTH_OF_GRAPH = 2      ///< DOG
  };
}

typedef chi_mesh::sweep_management::AngleSetGroup TAngleSetGroup;
typedef chi_mesh::sweep_management::AngleSet      TAngleSet;
typedef chi_mesh::sweep_management::STDG          TGSPO;
typedef std::vector<TGSPO>                        TLEVELED_GRAPH;

//###################################################################
class chi_mesh::sweep_management::SweepScheduler
{
private:
  SchedulingAlgorithm       scheduler_type;
  AngleAggregation&         angle_agg;

  struct RULE_VALUES
  {
    std::shared_ptr<TAngleSet> angle_set;
    int        depth_of_graph;
    int        sign_of_omegax;
    int        sign_of_omegay;
    int        sign_of_omegaz;
    size_t     set_index;

    explicit RULE_VALUES(std::shared_ptr<TAngleSet>& ref_as) :
      angle_set(ref_as)
    {
      depth_of_graph = 0;
      set_index      = 0;
      sign_of_omegax = 1;
      sign_of_omegay = 1;
      sign_of_omegaz = 1;
    }
  };
  std::vector<RULE_VALUES> rule_values;

  SweepChunk& m_sweep_chunk;
  const std::vector<size_t> sweep_timing_events_tag;

public:
  const size_t sweep_event_tag;

public:
  SweepScheduler(SchedulingAlgorithm in_scheduler_type,
                 AngleAggregation& in_angle_agg,
                 SweepChunk& in_sweep_chunk);

  void Sweep();
  double GetAverageSweepTime() const;
  std::vector<double> GetAngleSetTimings();
  SweepChunk& GetSweepChunk();

private:
  void ScheduleAlgoFIFO(SweepChunk& sweep_chunk);

  //02
  void InitializeAlgoDOG();
  void ScheduleAlgoDOG(SweepChunk& sweep_chunk);
};

#endif //CHI_SWEEPSCHEDULER_H