#include "geometry.h"


namespace sfm {

	/*Write the Coordinates*/
	void PointData::write(FILE *f)
	{
		fprintf(f, "%lf %lf %lf\n", pos[0], pos[1], pos[2]);
	}


	/* Read/write routines for tracks */
	void TrackData::read(FILE *f)
	{
		int size;
		fscanf(f, "%d", &size);

		for (int i = 0; i < size; i++) {
			ImageKey ik;
			fscanf(f, "%d %d", &(ik.first), &(ik.second));
			views.push_back(ik);
		}
	}

	void TrackData::write(FILE *f)
	{
		int size = (int)views.size();
		fprintf(f, "%d", size);

		for (int i = 0; i < size; i++) {
			fprintf(f, " %d %d", views[i].first, views[i].second);
		}

		fprintf(f, "\n");
	}

}