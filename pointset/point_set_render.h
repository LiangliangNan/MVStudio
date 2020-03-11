
#ifndef _POINT_SET_RENDERER_H_
#define _POINT_SET_RENDERER_H_



class PointSet;

class PointSetRender
{
public:
	PointSetRender();
	~PointSetRender(void);

	void set_pointset(PointSet* obj) { point_set_ = obj; }
	PointSet* target() const { return point_set_;  }

	float point_size() const  { return point_size_;  }
	void set_point_size(float x) { point_size_ = x;  }

	bool lighting() const { return lighting_; }
	void set_lighting(bool b) { lighting_ = b; }

	virtual void draw() ;

private:
	PointSet*	point_set_;
	float		point_size_;
	bool		lighting_;
};



#endif