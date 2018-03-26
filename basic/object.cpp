#include "object.h"
#include "canvas.h"
#include "logger.h"


Object::Object() {}


Object::~Object() {}


void Object::fit() {
	if (canvas_)
		canvas_->fit();
	else
		Logger::err("Object") << "you should assign canvas" << std::endl;
}


void Object::update_graphics() { 
	if (canvas_)	
		canvas_->update_graphics();
	else
		Logger::err("Object") << "you should assign canvas" << std::endl;
}

void Object::update_all() {
	if (canvas_)	
		canvas_->update_all();
	else
		Logger::err("Object") << "you should assign canvas" << std::endl;
}