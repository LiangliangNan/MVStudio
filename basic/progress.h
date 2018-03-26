#ifndef _BASIC_PROGRESS_H_
#define _BASIC_PROGRESS_H_

#include "basic_common.h"
#include "basic_types.h"


class ProgressClient;

/**
* For internal use, client code should use UserProgress.
*/
class BASIC_API Progress {
public:
	static void terminate();

	static Progress* instance();

	virtual void notify(std::size_t new_val, bool show_text = true);

	void set_client(ProgressClient* c) { client_ = c; }

	void push();
	void pop();

	void cancel() { canceled_ = true; }
	void clear_canceled() { canceled_ = false; }
	bool is_canceled() const { return canceled_; }

protected:
	Progress();
	virtual ~Progress();

	static Progress* instance_;
	ProgressClient* client_;
	int  level_;
	bool canceled_;
};

//_________________________________________________________

/**
* For internal use, client code do not need to use this one.
*/
class BASIC_API ProgressClient {
public:
	virtual void notify_progress(std::size_t new_val, bool show_text = true) = 0;
	virtual ~ProgressClient();
};

//_________________________________________________________

class BASIC_API ProgressLogger {
public:
	ProgressLogger(std::size_t max_val = 100, const std::string& task_name = "", bool quiet = false);
	virtual ~ProgressLogger();

	virtual void notify(std::size_t new_val, bool show_text = true);
	virtual void next();
	bool is_canceled() const {
		return Progress::instance()->is_canceled();
	}
	void reset(bool show_text = true) { notify(0, show_text); }
	void reset(std::size_t max_val, bool show_text = true);

protected:
	virtual void update(bool show_text = true);

private:
	std::size_t max_val_;
	std::string task_name_;
	std::size_t cur_val_;
	std::size_t cur_percent_;
	bool quiet_;
};


#endif

