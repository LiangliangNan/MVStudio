#include "camera.h"

#include "../basic/basic_types.h"
#include "../math/math_types.h"
#include "../math/matrix_driver.h"


namespace sfm {

	void get_intrinsics(const camera_params_t &camera, double *K) {
		if (!camera.known_intrinsics) {
			K[0] = camera.f;  K[1] = 0.0;       K[2] = 0.0;
			K[3] = 0.0;       K[4] = camera.f;  K[5] = 0.0;
			K[6] = 0.0;       K[7] = 0.0;       K[8] = 1.0;
		}
		else {
			memcpy(K, camera.K_known, 9 * sizeof(double));
		}
	}

	double camera_distance(camera_params_t* cam1, camera_params_t* cam2) {
		double center1[3];
		double Rinv1[9];
		matrix_invert(3, cam1->R, Rinv1);

		memcpy(center1, cam1->t, 3 * sizeof(double));

		double center2[3];
		double Rinv2[9];
		matrix_invert(3, cam2->R, Rinv2);

		memcpy(center2, cam2->t, 3 * sizeof(double));

		double dx = center1[0] - center2[0];
		double dy = center1[1] - center2[1];
		double dz = center1[2] - center2[2];

		return sqrt(dx * dx + dy * dy + dz * dz);
	}



	/* Compute the angle between two rays */
	double compute_ray_angle(vec2d p, vec2d q,
		const camera_params_t &cam1,
		const camera_params_t &cam2)
	{
		double K1[9], K2[9];
		get_intrinsics(cam1, K1);
		get_intrinsics(cam2, K2);

		double K1_inv[9], K2_inv[9];
		matrix_invert(3, K1, K1_inv);
		matrix_invert(3, K2, K2_inv);

		double p3[3] = { p.x, p.y, 1.0 };
		double q3[3] = { q.x, q.y, 1.0 };

		double p3_norm[3], q3_norm[3];
		matrix_product331(K1_inv, p3, p3_norm);
		matrix_product331(K2_inv, q3, q3_norm);

		vec2d p_norm(p3_norm[0] / p3_norm[2], p3_norm[1] / p3_norm[2]);
		vec2d q_norm(q3_norm[0] / q3_norm[2], q3_norm[1] / q3_norm[2]);

		double R1_inv[9], R2_inv[9];
		matrix_transpose(3, 3, (double *)cam1.R, R1_inv);
		matrix_transpose(3, 3, (double *)cam2.R, R2_inv);

		double p_w[3], q_w[3];

		double pv[3] = { p_norm.x, p_norm.y, -1.0 };
		double qv[3] = { q_norm.x, q_norm.y, -1.0 };

		double Rpv[3], Rqv[3];

		matrix_product331(R1_inv, pv, Rpv);
		matrix_product331(R2_inv, qv, Rqv);

		matrix_sum(3, 1, 3, 1, Rpv, (double *)cam1.t, p_w);
		matrix_sum(3, 1, 3, 1, Rqv, (double *)cam2.t, q_w);

		/* Subtract out the camera center */
		double p_vec[3], q_vec[3];
		matrix_diff(3, 1, 3, 1, p_w, (double *)cam1.t, p_vec);
		matrix_diff(3, 1, 3, 1, q_w, (double *)cam2.t, q_vec);

		/* Compute the angle between the rays */
		double dot;
		matrix_product(1, 3, 3, 1, p_vec, q_vec, &dot);

		double mag = matrix_norm(3, 1, p_vec) * matrix_norm(3, 1, q_vec);
		double ddm = dot / mag;
		ogf_clamp(ddm, -1.0 + 1.0e-8, 1.0 - 1.0e-8);
		return acos(ddm);
	}


	/* Check cheirality for a camera and a point */
	bool check_cheirality(vec3d p, const camera_params_t &camera)
	{
		double pt[3] = { p.x, p.y, p.z };
		double cam[3];

		pt[0] -= camera.t[0];
		pt[1] -= camera.t[1];
		pt[2] -= camera.t[2];
		matrix_product(3, 3, 3, 1, (double *)camera.R, pt, cam);

		// EDIT!!!
		if (cam[2] > 0.0)
			return false;
		else
			return true;
	}


	void invert_distortion(int n_in, int n_out, double r0, double r1,
		double *k_in, double *k_out)
	{
		const int NUM_SAMPLES = 20;

		int num_eqns = NUM_SAMPLES;
		int num_vars = n_out;

		double *A = new double[num_eqns * num_vars];
		double *b = new double[num_eqns];

		for (int i = 0; i < NUM_SAMPLES; i++) {
			double t = r0 + i * (r1 - r0) / (NUM_SAMPLES - 1);

			double p = 1.0;
			double a = 0.0;
			for (int j = 0; j < n_in; j++) {
				a += p * k_in[j];
				p = p * t;
			}

			double ap = 1.0;
			for (int j = 0; j < n_out; j++) {
				A[i * num_vars + j] = ap;
				ap = ap * a;
			}

			b[i] = t;
		}

		dgelsy_driver(A, b, k_out, num_eqns, num_vars, 1);

		delete[] A;
		delete[] b;
	}


	vec2d	 undistort_normalized_point(vec2d p, camera_params_t c) {
		double r = p.length();
		if (r == 0.0)
			return p;

		double t = 1.0;
		double a = 0.0;

		for (int i = 0; i < POLY_INVERSE_DEGREE; i++) {
			a += t * c.k_inv[i];
			t = t * r;
		}

		double factor = a / r;
		return (factor * p);

	}



	/* Finalize the camera */
	void Camera::Finalize() {
		/* Compute projection matrix */
		double K[9];
		double Ptmp[12] =
		{ R[0], R[1], R[2], t[0],
		R[3], R[4], R[5], t[1],
		R[6], R[7], R[8], t[2] };

		GetIntrinsics(K);

		matrix_product(3, 3, 3, 4, K, Ptmp, Pmatrix);
	}

	/* Get rigid transform for this camera */
	void Camera::GetRigid(double *T) const
	{
		T[0] = R[0];  T[1] = R[1];  T[2] = R[2];  T[3] = t[0];
		T[4] = R[3];  T[5] = R[4];  T[6] = R[5];  T[7] = t[1];
		T[8] = R[6];  T[9] = R[7];  T[10] = R[8]; T[11] = t[2];
	}

	/* Get a 4x4 rigid transform */
	void Camera::GetRigid4x4(double *T) const {
		GetRigid(T);
		T[12] = T[13] = T[14] = 0.0;
		T[15] = 1.0;
	}

	/* Set the position of the camera */
	void Camera::SetPosition(const double *pos) {
		matrix_product(3, 3, 3, 1, R, (double *)pos, t);
		matrix_scale(3, 1, t, -1.0, t);
	}

	/* Return the pose of the camera as a rotation matrix */
	void Camera::GetPose(double *R) const {
		matrix_transpose(3, 3, (double *)R, R);
	}

	void Camera::GetPoseQuaternion(double *q) const {
		double R[9];
		matrix_transpose(3, 3, (double *)R, R);
		matrix_to_quaternion(R, q);
	}

	/* Set the pose of the camera */
	void Camera::SetPose(const double *R) {
		matrix_transpose(3, 3, (double *)R, (double *)R);
	}

	/* Get upright rotation matrix */
	void Camera::GetUprightRotation(int rotate, double *R) {
		double R90[9] =
		{ 0.0, -1.0, 0.0,
		1.0, 0.0, 0.0,
		0.0, 0.0, 1.0 };

		double Rup[9];
		matrix_power(3, R90, rotate, Rup);

		matrix_product(3, 3, 3, 3, Rup, R, R);
	}

	/* Return the 3x3 intrinsic matrix */
	void Camera::GetIntrinsics(double *K) const {
		K[0] = focal; K[1] = 0.0;     K[2] = 0.0;
		K[3] = 0.0;     K[4] = focal; K[5] = 0.0;
		K[6] = 0.0;     K[7] = 0.0;     K[8] = 1.0;
	}

	/* Get the field of view */
	double Camera::GetFOV() const {
		return 2.0 * atan(width / (2.0 * focal));
	}

	double Camera::GetFOVMax(int rotate) const {
		if (((rotate % 2) == 0 && width >= height) ||
			((rotate % 2) == 1 && width < height)) {
			return 2.0 * atan(width / (2.0 * focal));
		}
		else {
			double vfov = 2.0 * atan(height / (2.0 * focal));
			double hfov = 2.0 * atan(tan(0.5 * vfov) * width / (double)height);
			printf("vfov = %0.3f, hfov = %0.3f\n", vfov, hfov);
			return hfov;
		}
	}

	/* Set the field of view */
	void Camera::SetFOV(double fov) {
		focal = 0.5 * width / tan(0.5 * DEG2RAD(fov));
	}

	/* Project the point into the camera */
	bool Camera::Project(const double *p, double *proj) const
	{
		double p4[4] = { p[0], p[1], p[2], 1.0 };
		double proj3[3];

		matrix_product(3, 4, 4, 1, (double *)Pmatrix, p4, proj3);

		if (proj3[2] == 0.0)
			return false;

		proj[0] = proj3[0] / -proj3[2];
		proj[1] = proj3[1] / -proj3[2];

		if (k[0] == 0.0 && k[1] == 0.0)
			return (proj3[2] < 0.0);

		double rsq = (proj[0] * proj[0] + proj[1] * proj[1]) / (focal * focal);
		double factor = 1.0 + k[0] * rsq + k[1] * rsq * rsq;

		if (rsq > 8.0 || factor < 0.0) // bad extrapolation
			return (proj[2] < 0.0);

		proj[0] *= factor;
		proj[1] *= factor;

		return (proj3[2] < 0.0);
	}

	/* Compute the essential matrix between this camera and another one */
	void Camera::ComputeEssentialMatrix(const Camera &cam,
		double *E, double *F)
	{
		/* Put the first camera at the canonical location */
		double P[16], Pinv[16];

		// memcpy(P, Pmatrix, sizeof(double) * 12);
		GetRigid(P);
		P[12] = P[13] = P[14] = 0.0;  P[15] = 1.0;

		matrix_invert(4, P, Pinv);

		double P2[16];
		cam.GetRigid(P2);
		P2[12] = P2[13] = P2[14] = 0.0;  P2[15] = 1.0;

		double P2new[16];
		matrix_product(4, 4, 4, 4, P2, Pinv, P2new);

		double R2new[9] = { P2new[0], P2new[1], P2new[2],
			P2new[4], P2new[5], P2new[6],
			P2new[8], P2new[9], P2new[10] };

		double t2new[3] = { P2new[3], P2new[7], P2new[11] };
		double t2new_cross[9] = { 0.0, -t2new[2], t2new[1],
			t2new[2], 0.0, -t2new[0],
			-t2new[1], t2new[0], 0.0 };

		matrix_product(3, 3, 3, 3, t2new_cross, R2new, E);

		/* Special black magic because we flipped the Z-axis */
		E[0] = -E[0];
		E[1] = -E[1];
		E[3] = -E[3];
		E[4] = -E[4];
		E[8] = -E[8];

		/* Compute F matrix */
		double K1[9], K2[9];
		GetIntrinsics(K1);
		cam.GetIntrinsics(K2);

		double K1inv[9], K2inv[9];
		matrix_invert(3, K1, K1inv);
		matrix_invert(3, K2, K2inv);

		double tmp[9];
		matrix_transpose_product(3, 3, 3, 3, K2inv, E, tmp);
		matrix_product(3, 3, 3, 3, tmp, K1inv, F);
	}

	/* Flip the camera over the z-axis */
	void Camera::Reflect()
	{
		R[2] = -R[2];
		R[5] = -R[5];
		R[6] = -R[6];
		R[7] = -R[7];

		t[2] = -t[2];

		Finalize();
	}

	/* Compute the distance to another camera */
	double Camera::CameraDistance(const Camera &cam) const
	{
		double pos1[3], pos2[3];

		GetPosition(pos1);
		cam.GetPosition(pos2);

		double diff[3];
		matrix_diff(3, 1, 3, 1, pos1, pos2, diff);

		return matrix_norm(3, 1, diff);
	}


	/* Returns true if the given 2D point is above the horizon */
	bool Camera::PointAboveHorizon(double *p)
	{
		double p3[3] = { p[0], p[1], 1.0 };
		double dot;
		matrix_product(1, 3, 3, 1, m_horizon, p3, &dot);

		return (dot > 0.0);
	}

	/* Returns true if the given point is in front of the camera */
	bool Camera::PointInFront(double *p)
	{
		/* Convert point to camera coordinates */
		double p_rot[3];
		matrix_product(3, 3, 3, 1, R, p, p_rot);
		matrix_sum(3, 1, 3, 1, p_rot, t, p_rot);

		/* Check z-coordinate */
		return (p_rot[2] < 0.0);
	}

	/* Interpolate between two camera views */
	Camera InterpolateCameras(const Camera &cam1,
		const Camera &cam2, double t)
	{
		double pos1[3], pos2[3];

		cam1.GetPosition(pos1);
		cam2.GetPosition(pos2);

		double pos_new[3];
		pos_new[0] = (1.0 - t) * pos1[0] + t * pos2[0];
		pos_new[1] = (1.0 - t) * pos1[1] + t * pos2[1];
		pos_new[2] = (1.0 - t) * pos1[2] + t * pos2[2];

		/* Get quaternions */
		double R1[9], R2[9];
		cam1.GetPose(R1);
		cam2.GetPose(R2);

		double q1[4], q2[4];
		matrix_to_quaternion(R1, q1);
		matrix_to_quaternion(R2, q2);

		double dot;
		matrix_product(1, 3, 3, 1, q1, q2, &dot);

		if (dot < 0.0) {
			matrix_scale(4, 1, q2, -1.0, q2);
		}

		double q_new[4];
		q_new[0] = (1.0 - t) * q1[0] + t * q2[0];
		q_new[1] = (1.0 - t) * q1[1] + t * q2[1];
		q_new[2] = (1.0 - t) * q1[2] + t * q2[2];
		q_new[3] = (1.0 - t) * q1[3] + t * q2[3];

		/* Normalize */
		double mag = matrix_norm(4, 1, q_new);
		matrix_scale(4, 1, q_new, 1.0 / mag, q_new);

		double Rnew[9];
		quaternion_to_matrix(q_new, Rnew);

		Camera cam_new;
		matrix_transpose(3, 3, Rnew, cam_new.R);

		matrix_product(3, 3, 3, 1, cam_new.R, pos_new, cam_new.t);
		matrix_scale(3, 1, cam_new.t, -1.0, cam_new.t);

		cam_new.focal = 1.0;

		return cam_new;
	}

	/* Convert a pixel position to a ray direction */
	void Camera::PixelToCameraRay(double x, double y, double *ray)
	{
		double ray_cam[3] = { x, y, -focal };

#if 0
		matrix_transpose_product(3, 3, 3, 1, R, ray_cam, ray);
#endif

		double norm = matrix_norm(3, 1, ray_cam);
		matrix_scale(3, 1, ray_cam, 1.0 / norm, ray);
	}

	/* Convert a pixel position to an (absolute) ray direction */
	void Camera::PixelToCameraRayAbsolute(double x, double y, double *ray)
	{
		double ray_cam[3] = { x, y, -focal };

		matrix_transpose_product(3, 3, 3, 1, R, ray_cam, ray);

		double norm = matrix_norm(3, 1, ray);
		matrix_scale(3, 1, ray, 1.0 / norm, ray);
	}

	static void SphToRot(double theta, double phi, double *R)
	{
		double v[3];

		v[0] = -cos(theta) * sin(phi);
		v[1] = -cos(phi);
		v[2] = -sin(theta) * sin(phi);

		/* Compute the new up vector */
		double phi_up = phi - 0.5 * M_PI;
		double theta_up = theta;
		double up[3];

		up[0] = cos(theta_up) * sin(phi_up);
		up[1] = cos(phi_up);
		up[2] = sin(theta_up) * sin(phi_up);

		double x_axis[3];
		matrix_cross(up, v, x_axis);

		memcpy(R + 0, x_axis, 3 * sizeof(double));
		memcpy(R + 3, up, 3 * sizeof(double));
		memcpy(R + 6, v, 3 * sizeof(double));
	}

	/* Point the camera in a different direction */
	void Camera::PointAt(double *ray)
	{
		double r = matrix_norm(3, 1, ray);

		/* Convert the ray to spherical coordinates */
		double theta = atan2(ray[2], ray[0]);
		double phi = acos(ray[1] / r);

		double R[9];
		SphToRot(theta, phi, R);

		double pos[3];
		GetPosition(pos);

		double Rnew[9];
		matrix_product(3, 3, 3, 3, R, R, Rnew);

		memcpy(R, Rnew, 9 * sizeof(double));
		SetPosition(pos);

		Finalize();
	}

	/* Point the camera in a different direction */
	void Camera::PointAtAbsolute(double *ray)
	{
		double r = matrix_norm(3, 1, ray);

		/* Convert the ray to spherical coordinates */
		double theta = atan2(ray[2], ray[0]);
		double phi = acos(ray[1] / r);

		double R[9];
		SphToRot(theta, phi, R);

		double pos[3];
		GetPosition(pos);

		memcpy(R, R, 9 * sizeof(double));
		SetPosition(pos);

		Finalize();
	}


	/* Return the view direction */
	void Camera::GetViewDirection(double *view) const
	{
		// double R[9];
		// GetPose(R);
		// double minuz[3] = { 0.0, 0.0, -1.0 };
		// matrix_product(3, 3, 3, 1, R, minuz, view);

		view[0] = -R[6];
		view[1] = -R[7];
		view[2] = -R[8];
	}

	/* Return the halfspace in front of the camera */
	void Camera::GetFrontHalfspace(double *plane) const
	{
		double v[3];
		GetViewDirection(v);

		double pos[3];
		GetPosition(pos);

		/* Set the position slightly in front of the real position */
		pos[0] += 1.0e-6 * v[0];
		pos[1] += 1.0e-6 * v[1];
		pos[2] += 1.0e-6 * v[2];

		double dot;
		matrix_product(1, 3, 3, 1, v, pos, &dot);

		plane[0] = v[0];
		plane[1] = v[1];
		plane[2] = v[2];
		plane[3] = -dot;
	}

	bool Camera::PointInsideImage(const double *p) const
	{
		double proj[2];
		bool in_front = Project(p, proj);
		if (!in_front) return false;
		return (proj[0] > -0.5*width && proj[0] < 0.5*width && proj[1]<0.5*height && proj[1]>-0.5*height);
	}

	/* Get a camera whose up vector points in the right direction */
	Camera Camera::GetUpCamera(double *up) const
	{
		double pos[3];
		GetPosition(pos);

		double up_image[3];
		matrix_product(3, 3, 3, 1, (double *)R, up, up_image);

		/* Project up-image onto the xy-plane */
		double up_image_proj[3] = { up_image[0], up_image[1], 0.0 };
		double mag = matrix_norm(3, 1, up_image_proj);
		matrix_scale(3, 1, up_image_proj, 1.0 / mag, up_image_proj);

		double yaxis[3] = { 0.0, 1.0, 0.0 };
		double dot;
		matrix_product(1, 3, 3, 1, up_image_proj, yaxis, &dot);

		double angle = acos(dot);

		double cross[3];
		matrix_cross(up_image_proj, yaxis, cross);

		mag = matrix_norm(3, 1, cross);
		matrix_scale(3, 1, cross, 1.0 / mag, cross);

		/* Rotate '-angle' radians around the cross axis */
		double Rroll[9];
		axis_angle_to_matrix(cross, -angle, Rroll);

		/* The world to camera transformation becomes Rroll * R */
#if 0
		double Rz[9] = { cos(angle), -sin(angle), 0.0,
			sin(angle),  cos(angle), 0.0,
			0.0,         0.0,        1.0   };

		double Rnew[9];
		matrix_product(3, 3, 3, 3, Rz, (double *) R, Rnew);
#endif

		double Rnew[9];
		matrix_transpose_product(3, 3, 3, 3, Rroll, (double *)R, Rnew);

		double RnewT[9];
		matrix_transpose(3, 3, Rnew, RnewT);

		Camera up_camera = *this;
		up_camera.SetPose(RnewT);
		up_camera.SetPosition(pos);
		up_camera.Finalize();

		matrix_product(3, 3, 3, 1, Rnew, up, up_image_proj);
		if (up_image_proj[0] > 1.0e-4)
			printf("Error in up camera computation\n");

		return up_camera;
	}

}