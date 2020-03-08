
#include "assertions.h"
#include "logger.h"

#include <cstdlib>


#ifndef WIN32
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#endif



static void ogf_abort() {
#ifdef WIN32
	std::abort() ;
#else
	int ogf_pid = getpid() ;
	Logger::err("Assert") << "Current pid is: " << ogf_pid << std::endl ;
	Logger::err("Assert") << "Going to bed" << std::endl ;
	Logger::err("Assert") << "Use: \'gdb MVStudio " << ogf_pid << "\' and then \'where\' to see the stack trace" << std::endl ;

	kill(ogf_pid, SIGSTOP) ;

	/*
	if(!fork()) {
	char spid[1024] ;
	sprintf("%d", spid, ogf_pid) ;
	execl("gdb", "gdb", "MVStudio", spid, nil) ;
	} else {
	kill(ogf_pid, SIGSTOP) ;
	}
	*/
#endif        
}

void ogf_assertion_failed(
						  const std::string& condition_string,
						  const std::string& file, int line
						  ) 
{
	Logger::err("Assert") << "Assertion failed: " << condition_string << std::endl ;
	Logger::err("Assert") << "File: " << file << std::endl ;
	Logger::err("Assert") << "Line: " << line << std::endl ;
	ogf_abort() ;
}

void ogf_range_assertion_failed(
								double value, double min_value, double max_value, 
								const std::string& file, int line
								) 
{
	Logger::err("Assert") << "Range assertion failed: " 
		<< value << " in " 
		<< "[ " << min_value << " ... " << max_value << " ]"
		<< std::endl ;
	Logger::err("Assert") << "File: " << file << std::endl ;
	Logger::err("Assert") << "Line: " << line << std::endl ;
	ogf_abort() ;
}

void ogf_should_not_have_reached(
								 const std::string& file, int line
								 ) 
{
	Logger::err("Assert") << "Control should not have reached this point:" << std::endl ;
	Logger::err("Assert") << "File: " << file << std::endl ;
	Logger::err("Assert") << "Line: " << line << std::endl ;
	ogf_abort() ;
}


