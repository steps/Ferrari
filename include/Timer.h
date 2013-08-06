//--------------------------------------------------------------------------------------------------
// Ferrari Reachability Index
// (c) 2012 Stephan Seufert. Web site: http://www.mpi-inf.mpg.de/~sseufert
//
// This work is licensed under the Creative Commons
// Attribution-Noncommercial-Share Alike 3.0 Unported License. To view a copy
// of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/
// or send a letter to Creative Commons, 171 Second Street, Suite 300,
// San Francisco, California, 94105, USA.
//--------------------------------------------------------------------------------------------------
#ifndef TIMER_H_
#define TIMER_H_
//--------------------------------------------------------------------------------------------------
#include <sys/time.h>
//--------------------------------------------------------------------------------------------------
class Timer {
private:
  struct timeval* _s;
public:
  Timer();
  ~Timer();
  void start();
  double stop();
};
//--------------------------------------------------------------------------------------------------
#endif /* TIMER_H_ */
