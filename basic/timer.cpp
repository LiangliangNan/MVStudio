
#include "timer.h"
#include "../basic/assertions.h"
#include "basic_types.h" // for "clip_precision()"

// Determine if the POSIX function getrusage is available, otherwise
// use the previous solution based on std::clock().
// First, detect POSIX. We cannot reliably use "unistd.h", 
// but limits.h is part of the C standard.
#include <climits>
#ifdef _POSIX_ARG_MAX // now that should be POSIX
#include <unistd.h>
#ifdef _POSIX_VERSION
#ifdef _XOPEN_UNIX // XSI: X/Open System Interfaces Extension
#define HAVE_GETRUSAGE 1
#endif
#endif
#endif


#ifdef HAVE_GETRUSAGE
// types, function prototype and constants for the POSIX function
// int getrusage (int who, struct rusage *usage);
#include <sys/resource.h>
// For the numerical limits
#else	//  HAVE_GETRUSAGE
// used for clock()
#include <ctime>
#endif	//  HAVE_GETRUSAGE

// For the numerical limits
#include <cfloat>



// Static member variable for Timer
// =====================================

bool Timer::failed_ = false;

// Member functions for Timer
// =====================================

double Timer::user_process_time() const {
	// Depends on the operating system.
	// Returns a (weakly ;-) monotone increasing time in seconds (with
	// possible wrap-around in case of overflow, see max()), or 0.0
	// if the system call for the time failed. If the system call
	// failed the static flag 'm_failed' is set and can be used
	// by the caller.
#ifdef HAVE_GETRUSAGE
	struct rusage usage;
	int ret = getrusage( RUSAGE_SELF, &usage);
	ogf_assert( ret == 0 );
	if ( ret == 0) {
		return double( usage.ru_utime.tv_sec)               // seconds
			+ double( usage.ru_utime.tv_usec) / 1000000.0;	// microseconds
	} else {
		//std::cerr << "Call to getrusage() in class Timer failed - timings will be 0." << std::endl;
        failed_ = true;
        return 0.0;
	}

#else // HAVE_GETRUSAGE
	std::clock_t clk = std::clock();
	ogf_assert( clk != (std::clock_t)-1 );
	if ( clk != (std::clock_t)-1 ) {
		return double(clk) / CLOCKS_PER_SEC;
	} else {
		//std::cerr << "Call to clock() in class Timer failed - timings will be 0." << std::endl;
        failed_ = true;
        return 0.0;
	}
#endif // HAVE_GETRUSAGE
}

double Timer::compute_precision() const {
	// Computes timer precision in seconds dynamically. Note that
	// the timer system call is probably non-trivial and will show
	// up in this time here (probably for one call). But that is just
	// fine that the call to the timer itself if reported as noise 
	// in the precision.
	double min_res = DBL_MAX;
	for ( int i = 0; i < 5; ++i) {
		double current = user_process_time();
		if ( failed_)
			return -1.0;
		double next    = user_process_time();
		while ( current >= next) { // wait until timer increases
			next = user_process_time();
			if ( failed_)
				return -1.0;
		}
		// Get the minimum timing difference of all runs.
		if ( min_res > next - current)
			min_res = next - current;
	}
	return min_res;
}

double Timer::precision() const {
	// computes precision upon first call
	// returns -1.0 if timer system call fails.
	static double prec = compute_precision();
	return prec;
}


// -----------------------------------------------------------------------

// Member functions for Timer
// ===========================

void Timer::start() {
	ogf_assert( ! running_);
	started_ = user_process_time();
	running_ = true;
	++ interv_;
}

void Timer::stop() {
	ogf_assert( running_);
	double t = user_process_time();
	elapsed_ += (t - started_);
	started_ = 0.0;
	running_ = false;
}

void Timer::reset() {
	interv_  = 0;
	elapsed_ = 0.0;
	if (running_) {
		started_ = user_process_time();
		++ interv_;
	} else {
		started_ = 0.0;
	}
}

double Timer::time() const {
	if (running_) {
		double t = user_process_time();
		return clip_precision(elapsed_ + (t - started_), 2);
		//return elapsed_ + (t - started_);
	}

	return clip_precision(elapsed_, 2);
	//return elapsed_;
}
