
#ifndef _STOP_WATCH_H_
#define _STOP_WATCH_H_

#include "basic_common.h"


#ifdef WIN32
#	include <windows.h>
#else 
#	include <sys/time.h>
#endif // WIN32

//______________________________________________________________________


/**
* 
* The purpose of this file is to make a timer function
* that is as precise as possible on any given platform.
*
* usage example:
*   {
*      StopWatch w ;
*      // do task_1 ...
*      std::cout << "task_1 done. time: " << w.elapsed() << " seconds";
*	   w.start();
*      // do task_2 ...
*      std::cout << "task_2 done. time: " << w.elapsed() << " seconds";
*   } 
*/

class BASIC_API StopWatch 
{
public :
	StopWatch() ; // the watch will automatically start in construction
	~StopWatch() ;

	void  start() ;

	// returns user elapsed time since the construction / start in seconds.
	double elapsed() const ;

private:

#ifdef WIN32
	LONGLONG  freq_;
	LONGLONG  start_count_;
#else
	timeval start_time_;
#endif

} ;


#endif

