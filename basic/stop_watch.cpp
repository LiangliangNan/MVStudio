
#include "stop_watch.h"
#include "basic_types.h" // for "clip_precision()"
#include <iostream>



//_________________________________________________________

StopWatch::StopWatch() {
	start();
}

StopWatch::~StopWatch() {}


void StopWatch::start() {
#ifdef WIN32
	LARGE_INTEGER  largeInteger;
	QueryPerformanceFrequency(&largeInteger);
	freq_ = largeInteger.QuadPart;
	QueryPerformanceCounter(&largeInteger);
	start_count_ = largeInteger.QuadPart;
#else
	gettimeofday(&start_time_, 0);
#endif // WIN32
}

double StopWatch::elapsed() const {
#ifdef WIN32
	LARGE_INTEGER  largeInteger;
	QueryPerformanceCounter(&largeInteger);
	LONGLONG now_count = largeInteger.QuadPart;
	double time = (double)( (now_count - start_count_) / static_cast<double>(freq_) );
	return clip_precision(time, 2);
#else
	timeval now;
	gettimeofday(&now, 0);
    double time = (now.tv_sec - start_time_.tv_sec) + (now.tv_usec - start_time_.tv_usec) / 1.0e6;
    return clip_precision(time, 2);
#endif  // WIN32
}
