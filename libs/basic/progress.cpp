
#include "progress.h"
#include "assertions.h"


Progress* Progress::instance_ = nil;


Progress* Progress::instance() {
	if (instance_ == nil)
		instance_ = new Progress;
	return instance_;
}

void Progress::terminate() {
	//mpl_debug_assert(instance_ != nil);
	if (instance_) {
		delete instance_;
		instance_ = nil;
	}
}

Progress::Progress() : client_(nil), level_(0), canceled_(false) {
}

Progress::~Progress() {
}

void Progress::push() {
	level_++;
	if (level_ == 1) {
		clear_canceled();
	}
}

void Progress::pop() {
	ogf_debug_assert(level_ > 0);
	level_--;
}

void Progress::notify(std::size_t new_val, bool show_text) {
	if (client_ != nil && level_ < 2) {
		client_->notify_progress(new_val, show_text);
	}
}

//_________________________________________________________


ProgressClient::~ProgressClient() {
	// the Progress will take care of the destroying itself
	// 	if (Progress::instance())
	// 		delete Progress::instance();
}

//_________________________________________________________

ProgressLogger::ProgressLogger(std::size_t max_val, const std::string& task_name, bool quiet)
	: max_val_(max_val), task_name_(task_name), quiet_(quiet)
{
	cur_val_ = 0;
	cur_percent_ = 0;
	Progress::instance()->push();
	if (!quiet_) {
		Progress::instance()->notify(0, true);
	}
}

void ProgressLogger::reset(std::size_t max_val, bool show_text) {
	max_val_ = max_val;
	reset(show_text);
}

ProgressLogger::~ProgressLogger() {
	Progress::instance()->notify(100, true);
	Progress::instance()->notify(0, false);
	Progress::instance()->pop();
}

void ProgressLogger::next() {
	cur_val_++;
	update(true);
}

void ProgressLogger::notify(std::size_t new_val, bool show_text) {
	cur_val_ = new_val;
	update(show_text);
}


void ProgressLogger::update(bool show_text) {
	std::size_t percent = cur_val_ * 100 / ogf_max<std::size_t>(1, max_val_ - 1);
	if (percent != cur_percent_) {
		cur_percent_ = percent;
		if (!quiet_) {
			Progress::instance()->notify(ogf_min<std::size_t>(cur_percent_, 100), show_text);
		}
	}
}

