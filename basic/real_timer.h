
#ifndef _BASIC_REAL_TIMER_H_
#define _BASIC_REAL_TIMER_H_

#include "basic_common.h"
#include <cfloat>	// For the numerical limits


// Liangliang: this implementation is modified from CGAL. 
// See <CGAL/Real_timer.h> and <CGAL/Real_timer.cpp>.

// SECTION: A Timer Measuring Real-Time
// ========================================================================
// 
// DEFINITION
// 
// The class RealTimer is a timer class for measuring real time. A timer 't' 
// of type RealTimer is an object with a state. It is either running or it 
// is stopped. The state is controlled with t.start() and t.stop(). The timer 
// counts the time elapsed since its creation or last reset. It counts only 
// the time where it is in the running state. The time information is given 
// in seconds. The timer counts also the number of intervals it was running, 
// i.e. it counts the number of calls of the start() member function since 
// the last reset. If the reset occurs while the timer is running it counts 
// as the first interval

class BASIC_API RealTimer 
{
public:
	RealTimer() : elapsed_(0.0), started_(0.0), interv_(0), running_(false) {}

	void	start();
	void	stop ();
	void	reset();
	bool	is_running() const { return running_; }

	double	time()       const;
	int		intervals()  const { return interv_; }

	// Returns timer precision. Computes it dynamically at first call.
	// Returns -1.0 if timer system call fails, which, for a proper coded
	// test towards precision leads to an immediate stop of an otherwise 
	// infinite loop (fixed tolerance * total time >= precision).
	double	precision()  const;

private:
	double	get_real_time()     const; // in seconds
	double	compute_precision() const; // in seconds

	double	elapsed_;
	double	started_;
	int		interv_;
	bool	running_;

	static bool failed_;
};





#endif