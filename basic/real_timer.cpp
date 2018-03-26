
#include "real_timer.h"
#include "../basic/assertions.h"
#include "basic_types.h" // for "clip_precision()"


#if defined (_MSC_VER)
#  include <sys/timeb.h>
#  include <sys/types.h>

#elif defined (__MINGW32__)
#  include <sys/timeb.h>
#  include <sys/types.h>

#else
// If none of the above PC compilers, use POSIX fct. gettimeofday()
#  include <sys/time.h>
#endif



// Static member variable for RealTimer
// =====================================

bool RealTimer::failed_ = false;

// Member functions for RealTimer
// =====================================

double RealTimer::get_real_time() const {
	// Depends on the operating system.
	// Returns a (weakly ;-) monotone increasing time in seconds (with
	// possible wrap-around in case of overflow, see max()), or 0.0
	// if the system call for the time failed. If the system call
	// failed the static flag 'failed_' is set and can be used
	// by the caller.
#if   defined(_MSC_VER)
	struct _timeb  t;
	_ftime_s(&t);  
	return double(t.time) + double(t.millitm) / 1000.0;
#elif defined (__MINGW32__)
	struct timeb t;
	ftime(&t);
	return double(t.time) + double(t.millitm) / 1000.0;
#else // ! _MSC_VER && ! __MINGW32__//
	struct timeval t;
	int ret = gettimeofday( &t, NULL);
	ogf_assert(ret == 0);
	if ( ret == 0) {
		return double(t.tv_sec) + double(t.tv_usec) / 1000000;
	} else {
		//std::cerr << "Call to gettimeofday() in class RealTimer failed - timings will be 0." << std::endl;
        failed_ = true;
        return 0.0;
	}
#endif // ! _MSC_VER && ! __MINGW32__//
}

double RealTimer::compute_precision() const {
	// Computes timer precision in seconds dynamically. Note that
	// the timer system call is probably non-trivial and will show
	// up in this time here (probably for one call). But that is just
	// fine that the call to the timer itself if reported as noise 
	// in the precision.
	double min_res = DBL_MAX;
	for ( int i = 0; i < 5; ++i) {
		double current = get_real_time();
		if ( failed_)
			return -1.0;
		double next    = get_real_time();
		while ( current >= next) { // wait until timer increases
			next = get_real_time();
			if ( failed_)
				return -1.0;
		}
		// Get the minimum timing difference of all runs.
		if ( min_res > next - current)
			min_res = next - current;
	}
	return min_res;
}

double RealTimer::precision() const {
	// computes precision upon first call
	// returns -1.0 if timer system call fails.
	static double prec = compute_precision();
	return prec;
}


// -----------------------------------------------------------------------

// Member functions for RealTimer
// ===========================

void RealTimer::start() {
	ogf_assert( ! running_);
	started_ = get_real_time();
	running_ = true;
	++ interv_;
}

void RealTimer::stop() {
	ogf_assert( running_);
	double t = get_real_time();
	elapsed_ += (t - started_);
	started_ = 0.0;
	running_ = false;
}

void RealTimer::reset() {
	interv_  = 0;
	elapsed_ = 0.0;
	if (running_) {
		started_ = get_real_time();
		++ interv_;
	} else {
		started_ = 0.0;
	}
}

double RealTimer::time() const {
	if (running_) {
		double t = get_real_time();
		return clip_precision(elapsed_ + (t - started_), 2);
		//return elapsed_ + (t - started_);
	}

	return clip_precision(elapsed_, 2);
	//return elapsed_;
}
