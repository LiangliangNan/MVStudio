
#ifndef _SFM_CAMERA_H_
#define _SFM_CAMERA_H_

#include <easy3d/core/types.h>

using namespace easy3d;

namespace sfm {

#define NUM_CAMERA_PARAMS 9
#define POLY_INVERSE_DEGREE 6

	typedef struct {
		double R[9];     /* Rotation */
		double t[3];     /* Translation */
		double f;        /* Focal length */
		double k[2];     /* Undistortion parameters */
		double k_inv[POLY_INVERSE_DEGREE]; /* Inverse undistortion parameters */
		char constrained[NUM_CAMERA_PARAMS];
		double constraints[NUM_CAMERA_PARAMS];  /* Constraints (if used) */
		double weights[NUM_CAMERA_PARAMS];      /* Weights on the constraints */
		double K_known[9];  /* Intrinsics (if known) */
		double k_known[5];  /* Distortion params (if known) */

		char known_intrinsics;   /* Are the intrinsics known? */

		double f_scale, k_scale; /* Scale on focal length, distortion params */
	} camera_params_t;


	void	 get_intrinsics(const camera_params_t &camera, double *K);

	double camera_distance(const camera_params_t* cam1, const camera_params_t* cam2);

	/* Compute the angle between two rays */
	double compute_ray_angle(easy3d::dvec2 p, easy3d::dvec2 q, const camera_params_t &cam1, const camera_params_t &cam2);

	/* Check cheirality for a camera and a point */
	bool	 check_cheirality(easy3d::dvec3 p, const camera_params_t &camera);

	easy3d::dvec2	 undistort_normalized_point(easy3d::dvec2 p, camera_params_t c);

	void   invert_distortion(int n_in, int n_out, double r0, double r1,
		double *k_in, double *k_out);






	class Camera {
	public:
		Camera() {
			adjusted = false;

			constrained[0] = constrained[1] = constrained[2] = false;
			constrained[3] = constrained[4] = constrained[5] = false;
			constrained[6] = false;

			constraints[0] = constraints[1] = constraints[2] = 0.0;
			constraints[3] = constraints[4] = constraints[5] = 0.0;
			constraints[6] = 0.0;

			constraint_weights[0] = constraint_weights[1] = constraint_weights[2] = 0.0;
			constraint_weights[3] = constraint_weights[4] = constraint_weights[5] = 0.0;
			constraint_weights[6] = 0.0;

			k[0] = k[1] = 0.0;
		}

		/* Finalize the camera */
		void Finalize();
		/* Get rigid transform for this camera */
		void GetRigid(double *T) const;
		/* Get a 4x4 rigid transform */
		void GetRigid4x4(double *T) const;
		/* Return the position of the camera */
		void inline GetPosition(double *pos) const {
			pos[0] = R[0] * t[0] + R[3] * t[1] + R[6] * t[2];
			pos[1] = R[1] * t[0] + R[4] * t[1] + R[7] * t[2];
			pos[2] = R[2] * t[0] + R[5] * t[1] + R[8] * t[2];

			pos[0] = -pos[0];
			pos[1] = -pos[1];
			pos[2] = -pos[2];
		}

		/* Set the position of the camera */
		void SetPosition(const double *pos);
		/* Return the pose of the camera as a rotation matrix */
		void GetPose(double *R) const;
		void GetPoseQuaternion(double *q) const;
		/* Set the pose of the camera */
		void SetPose(const double *R);
		/* Get upright rotation matrix */
		void GetUprightRotation(int rotate, double *R);
		/* Return the 3x3 intrinsic matrix */
		void GetIntrinsics(double *K) const;
		/* Get the field of view */
		double GetFOV() const;
		double GetFOVMax(int rotate) const;
		/* Set the field of view */
		void SetFOV(double fov);
		/* Project the point into the camera */
		bool Project(const double *p, double *proj) const;
		/* Compute the essential matrix between this camera and another one */
		void ComputeEssentialMatrix(const Camera &cam, double *E, double *F);
		/* Flip the camera over the z-axis */
		void Reflect();
		/* Compute the distance to another camera */
		double CameraDistance(const Camera &cam) const;
		/* Returns true if the given point is in front of the camera */
		bool PointInFront(double *p);
		/* Returns true if the given 2D point is above the horizon */
		bool PointAboveHorizon(double *p);

		/* Convert a pixel position to a ray direction */
		void PixelToCameraRay(double x, double y, double *ray);
		/* Convert a pixel position to an (absolute) ray direction */
		void PixelToCameraRayAbsolute(double x, double y, double *ray);
		/* Point the camera in a different direction */
		void PointAt(double *ray);
		void PointAtAbsolute(double *ray);
		/* Return the view direction */
		void GetViewDirection(double *view) const;
		/* Return the halfspace in front of the camera */
		void GetFrontHalfspace(double *plane) const;
		bool PointInsideImage(const double *p) const;

		/* Get a camera whose up vector points in the right direction */
		Camera GetUpCamera(double *up) const;

		bool		adjusted;		/* Has this camera been adjusted? */
		double	focal;		/* Focal length */
		double	k[2];			/* Distortion parameters */
		double	R[9], t[3];		/* Extrinsics */
		double	Pmatrix[12];
		int		width, height;	/* Image dimensions */

		/* Horizon line */
		double m_horizon[3];

		/* Constraints on camera center location */
		bool		constrained[7];
		double	constraints[7];
		double	constraint_weights[7];

	};

	/* Interpolate between two camera views */
	Camera interpolate_cameras(const Camera &cam1, const Camera &cam2, double t);


}


#endif /* __camera_h__ */
