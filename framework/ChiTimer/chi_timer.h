#ifndef CHI_TIMER_H
#define CHI_TIMER_H

#include <string>
#include <chrono>


//################################################################### CLASS DEF
namespace chi_objects
{
  /** Timer object.*/
  class ChiTimer
  {
  public:
    std::chrono::steady_clock::time_point startTime;

  public:
    //00
                ChiTimer() noexcept;
    //01
    void   		  Reset();
    double 		  GetTime() const;
    std::string GetTimeString() const;
    static std::string GetLocalDateTimeString();
  };
}

#endif

