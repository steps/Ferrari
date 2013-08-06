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
#include "Timer.h"
#include <cassert>
//--------------------------------------------------------------------------------------------------
Timer::Timer() : _s(0) {
}
//--------------------------------------------------------------------------------------------------
Timer::~Timer() {
  if (_s) delete _s;
}
//--------------------------------------------------------------------------------------------------
void Timer::start() {
  assert(!_s);
  _s = new timeval();
  gettimeofday(_s, 0);
}
//--------------------------------------------------------------------------------------------------
double Timer::stop() {
  assert(_s);
  struct timeval* e = new timeval();
  gettimeofday(e, 0);
  long ms_s = _s->tv_sec * 1000000 + _s->tv_usec;
  long ms_e = e->tv_sec * 1000000 + e->tv_usec;
  double time = (ms_e - ms_s) / 1000.0;
  delete e;
  delete _s;
  _s = 0;
  return time;
}
//--------------------------------------------------------------------------------------------------

