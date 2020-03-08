
#include "5point.h"
#include "poly3.h"
#include "fmatrix.h"
#include "triangulate.h"

#include "../math/matrix_driver.h"

#include <algorithm>


using namespace easy3d;

void compute_nullspace_basis(int n, dvec2* a, dvec2* b, double *basis) {
	if (n < 5) {
		fprintf(stderr, "[compute_nullspace_basis] n must be >= 5\n");
		return;
	}

	int max_dim = std::max(9, n);

	double* Q = new double[max_dim * max_dim];
	double* S = new double[max_dim];
	double* U = new double[max_dim * max_dim];

	/* Create the 5x9 epipolar constraint matrix */
	for (int i = 0; i < 5; i++) {
		double *row = Q + i * 9;

		row[0] = a[i].x * b[i].x;
		row[1] = a[i].y * b[i].x;
		row[2] = b[i].x;

		row[3] = a[i].x * b[i].y;
		row[4] = a[i].y * b[i].y;
		row[5] = b[i].y;

		row[6] = a[i].x;
		row[7] = a[i].y;
		row[8] = 1.0;
	}

	/* Find four vectors that span the right nullspace of the matrix */
	double VT[81];
	dgesvd_driver(n, 9, Q, U, S, VT);

	memcpy(basis, VT + 5 * 9, 36 * sizeof(double));

	delete[] Q;
	delete[] S;
	delete[] U;
}

void compute_constraint_matrix(double *basis, poly3_t *constraints)
{
    /* Basis rows are X, Y, Z, W 
     * Essential matrix is or form x*X + y*Y + z*Z + W */

    /* Create a polynomial for each entry of E */
    poly3_t polys[9];
    poly3_t poly_term1, poly_term2, poly_term3, poly_det;
    poly3_t poly_EET[6], poly_lambda[6], poly_tr, poly_lambdaE[9];

    int i;

    for (i=0; i < 9; i++) {
        polys[i] = poly3_new(basis[i], basis[9+i], basis[18+i], basis[27+i]);
    }

    /* Create a polynormial from the constraint det(E) = 0 */
    poly_term1 = poly3_mult21( poly3_sub( poly3_mult11(polys[1], polys[5]), 
                                          poly3_mult11(polys[2], polys[4]) ), 
                               polys[6] );

    // matrix_print(1, 20, poly_term1.v);
    
    poly_term2 = poly3_mult21( poly3_sub( poly3_mult11(polys[2], polys[3]), 
                                          poly3_mult11(polys[0], polys[5]) ), 
                               polys[7] );

    poly_term3 = poly3_mult21( poly3_sub( poly3_mult11(polys[0], polys[4]), 
                                          poly3_mult11(polys[1], polys[3]) ), 
                               polys[8] );

    poly_det = poly3_add(poly_term1, poly3_add(poly_term2, poly_term3));
    
    /* Create polynomials for the singular value constraint */
    for (i = 0; i < 6; i++) {
        int r = 0, c = 0, k;
        poly_EET[i] = poly3_new(0.0, 0.0, 0.0, 0.0);
        switch(i) {
        case 0:
        case 1:
        case 2:
            r = 0;
            c = i;
            break;
        case 3:
        case 4:
            r = 1;
            c = i-2;
            break;
        case 5:
            r = 2;
            c = 2;
            break;
        }

        for (k = 0; k < 3; k++) {
            poly_EET[i] = poly3_add(poly_EET[i], poly3_mult11(polys[r*3+k], polys[c*3+k]));
        }
    }

    poly_tr = poly3_add3(poly_EET[0], poly_EET[3], poly_EET[5]);
    poly_tr = poly3_scale(poly_tr, 0.5);

    poly_lambda[0] = poly3_sub(poly_EET[0], poly_tr);
    poly_lambda[1] = poly_EET[1];
    poly_lambda[2] = poly_EET[2];
    
    poly_lambda[3] = poly3_sub(poly_EET[3], poly_tr);
    poly_lambda[4] = poly_EET[4];
    
    poly_lambda[5] = poly3_sub(poly_EET[5], poly_tr);
    
    poly_lambdaE[0] = poly3_add3(poly3_mult(poly_lambda[0], polys[0]),
                                 poly3_mult(poly_lambda[1], polys[3]),
                                 poly3_mult(poly_lambda[2], polys[6]));
    
    poly_lambdaE[1] = poly3_add3(poly3_mult(poly_lambda[0], polys[1]),
                                 poly3_mult(poly_lambda[1], polys[4]),
                                 poly3_mult(poly_lambda[2], polys[7]));

    poly_lambdaE[2] = poly3_add3(poly3_mult(poly_lambda[0], polys[2]),
                                 poly3_mult(poly_lambda[1], polys[5]),
                                 poly3_mult(poly_lambda[2], polys[8]));

    poly_lambdaE[3] = poly3_add3(poly3_mult21(poly_lambda[1], polys[0]),
                                 poly3_mult21(poly_lambda[3], polys[3]),
                                 poly3_mult21(poly_lambda[4], polys[6]));
    
    poly_lambdaE[4] = poly3_add3(poly3_mult21(poly_lambda[1], polys[1]),
                                 poly3_mult21(poly_lambda[3], polys[4]),
                                 poly3_mult21(poly_lambda[4], polys[7]));

    poly_lambdaE[5] = poly3_add3(poly3_mult21(poly_lambda[1], polys[2]),
                                 poly3_mult21(poly_lambda[3], polys[5]),
                                 poly3_mult21(poly_lambda[4], polys[8]));

    poly_lambdaE[6] = poly3_add3(poly3_mult21(poly_lambda[2], polys[0]),
                                 poly3_mult21(poly_lambda[4], polys[3]),
                                 poly3_mult21(poly_lambda[5], polys[6]));
    
    poly_lambdaE[7] = poly3_add3(poly3_mult21(poly_lambda[2], polys[1]),
                                 poly3_mult21(poly_lambda[4], polys[4]),
                                 poly3_mult21(poly_lambda[5], polys[7]));

    poly_lambdaE[8] = poly3_add3(poly3_mult21(poly_lambda[2], polys[2]),
                                 poly3_mult21(poly_lambda[4], polys[5]),
                                 poly3_mult21(poly_lambda[5], polys[8]));
 
    for (i=0; i < 9; i++) 
        constraints[i] = poly_lambdaE[i];

    constraints[9] = poly_det;
}

void compute_Grabner_basis(poly3_t *constraints, double *Gbasis) 
{
    double A[200];
    int i, j;

    for (i = 0; i < 10; i++) {
        // memcpy(A + 20 * i, constraints[i].v, sizeof(double) * 20);
        double *row = A + 20 * i;

        /* x3 x2y xy2 y3 x2z xyz y2z xz2 yz2 z3 x2 xy y2 xz yz z2 x  y  z  1 */

        row[0] = constraints[i].v[POLY3_X3];
        row[1] = constraints[i].v[POLY3_X2Y];
        row[2] = constraints[i].v[POLY3_XY2];
        row[3] = constraints[i].v[POLY3_Y3];
        row[4] = constraints[i].v[POLY3_X2Z];
        row[5] = constraints[i].v[POLY3_XYZ];
        row[6] = constraints[i].v[POLY3_Y2Z];
        row[7] = constraints[i].v[POLY3_XZ2];
        row[8] = constraints[i].v[POLY3_YZ2];
        row[9] = constraints[i].v[POLY3_Z3];
        row[10] = constraints[i].v[POLY3_X2];
        row[11] = constraints[i].v[POLY3_XY];
        row[12] = constraints[i].v[POLY3_Y2];
        row[13] = constraints[i].v[POLY3_XZ];
        row[14] = constraints[i].v[POLY3_YZ];
        row[15] = constraints[i].v[POLY3_Z2];
        row[16] = constraints[i].v[POLY3_X];
        row[17] = constraints[i].v[POLY3_Y];
        row[18] = constraints[i].v[POLY3_Z];
        row[19] = constraints[i].v[POLY3_UNIT];
    }
    
    /* Do a full Gaussian elimination */

    for (i = 0; i < 10; i++) {
        /* Make the leading coefficient of row i = 1 */
        double leading = A[20 * i + i];
        matrix_scale(20, 1, A + 20 * i, 1.0 / leading, A + 20 * i);

        /* Subtract from other rows */
        for (j = i+1; j < 10; j++) {
            double leading2 = A[20 * j + i];
            double scaled_row[20];
            matrix_scale(20, 1, A + 20 * i, leading2, scaled_row);
            matrix_diff(20, 1, 20, 1, A + 20 * j, scaled_row, A + 20 * j);
        }
    }

    /* Now, do the back substitution */
    for (i = 9; i >= 0; i--) {
        for (j = 0; j < i; j++) {
            double scale = A[20 * j + i];
            double scaled_row[20];
            matrix_scale(20, 1, A + 20 * i, scale, scaled_row);
            matrix_diff(20, 1, 20, 1, A + 20 * j, scaled_row, A + 20 * j);
        }
    }
    
    /* Copy out results */
    for (i = 0; i < 10; i++) {
        memcpy(Gbasis + i * 10, A + i * 20 + 10, sizeof(double) * 10);
    }
}

void compute_action_matrix(double *Gbasis, double *At) 
{
    int i;
    for (i = 0; i < 100; i++)
        At[i] = 0.0;
    
    matrix_scale(10, 1, Gbasis +  0, -1.0, At +  0);
    matrix_scale(10, 1, Gbasis + 10, -1.0, At + 10);
    matrix_scale(10, 1, Gbasis + 20, -1.0, At + 20);
    matrix_scale(10, 1, Gbasis + 40, -1.0, At + 30);
    matrix_scale(10, 1, Gbasis + 50, -1.0, At + 40);
    matrix_scale(10, 1, Gbasis + 70, -1.0, At + 50);

    At[6 * 10 + 0] = 1.0;
    At[7 * 10 + 1] = 1.0;
    At[8 * 10 + 3] = 1.0;
    At[9 * 10 + 6] = 1.0;
}

void compute_Ematrices_Gb(double *At, double *basis, int *num_solns, double *E)
{
    double evec[100], eval[10];
    int real = dgeev_driver(10, At, evec, eval);
    int i;
    
    for (i = 0; i < real; i++) {
        double X[9], Y[9], Z[9], tmp1[9], tmp2[9];

        double x = evec[10 * i + 6];
        double y = evec[10 * i + 7];
        double z = evec[10 * i + 8];
        double w = evec[10 * i + 9];
        double w_inv = 1.0 / w;
        
        x = x * w_inv;
        y = y * w_inv;
        z = z * w_inv;

        matrix_scale(9, 1, basis + 0, x, X);
        matrix_scale(9, 1, basis + 9, y, Y);
        matrix_scale(9, 1, basis + 18, z, Z);

        matrix_sum(9, 1, 9, 1, X, Y, tmp1);
        matrix_sum(9, 1, 9, 1, Z, basis + 27, tmp2);
        matrix_sum(9, 1, 9, 1, tmp1, tmp2, E + 9 * i);

        matrix_scale(9, 1, E + 9 * i, 1.0 / E[9 * i + 8], E + 9 * i);
    }

    *num_solns = real;
}

void generate_Ematrix_hypotheses(int n, dvec2* rt_pts, dvec2* left_pts,
                                 int *num_poses, double *E) 
{
    double basis[36], Gbasis[100], At[100];
    poly3_t constraints[10];

    if (n < 5) {
        fprintf(stderr, "[generate_Ematrix_hypotheses] n must be >= 5\n");
        return;
    }

    /* Generate the nullspace basis of the epipolar constraint matrix */
    compute_nullspace_basis(n, rt_pts, left_pts, basis);
    compute_constraint_matrix(basis, constraints);

    compute_Grabner_basis(constraints, Gbasis);
    compute_action_matrix(Gbasis, At);
    compute_Ematrices_Gb(At, basis, num_poses, E);
}

void choose(int n, int k, int *arr)
{
    int i;
    
    if (k > n) {
        fprintf(stderr, "[choose] Error: k > n\n");
        return;
    }

    for (i = 0; i < k; i++) {
        while (1) {
            int idx = rand() % n;
            int j, redo = 0;

            for (j = 0; j < i; j++) {
                if (idx == arr[j]) {
                    redo = 1;
                    break;
                }
            }

            if (!redo) {
                arr[i] = idx;
                break;
            }
        }
    }
}

int evaluate_Ematrix(int n, dvec2* r_pts, dvec2* l_pts, double thresh_norm,
                     double *F, int *best_inlier, double *score)
{
    int num_inliers = 0;
    int i;
    double min_resid = 1.0e20;
    double likelihood = 0.0;

    for (i = 0; i < n; i++) {
        dvec3 r(r_pts[i].x, r_pts[i].y, 1.0);
		dvec3 l(l_pts[i].x, l_pts[i].y, 1.0);

        double resid = fmatrix_compute_residual(F, l, r);    
   
        likelihood += log(1.0 + resid * resid / (thresh_norm));

        if (resid < thresh_norm) {
            num_inliers++;

            if (resid < min_resid) {
                min_resid = resid;
                *best_inlier = i;
            }
        }
    }

    *score = likelihood;
    // *score = 1.0 / num_inliers;

    return num_inliers;
}

int compute_pose_ransac(int n, dvec2* r_pts, dvec2* l_pts,
                        double *K1, double *K2, 
                        double ransac_threshold, int ransac_rounds, 
                        double *R_out, double *t_out)
{
    int i, round;
    double thresh_norm;
    double K1_inv[9], K2_inv[9];
    int max_inliers = 0;
    double min_score = DBL_MAX;
    double E_best[9];
    dvec2 r_best, l_best;

	dvec2* r_pts_norm = new dvec2[n];
	dvec2* l_pts_norm = new dvec2[n];

    matrix_invert(3, K1, K1_inv);
    matrix_invert(3, K2, K2_inv);

    for (i = 0; i < n; i++) {
        double r[3] = { r_pts[i].x, r_pts[i].y, 1.0 };
        double l[3] = { l_pts[i].x, l_pts[i].y, 1.0 };

        double r_norm[3], l_norm[3];

        matrix_product331(K1_inv, r, r_norm);
        matrix_product331(K2_inv, l, l_norm);

        r_pts_norm[i] = dvec2(-r_norm[0], -r_norm[1]);
		l_pts_norm[i] = dvec2(-l_norm[0], -l_norm[1]);
    }

    thresh_norm = ransac_threshold * ransac_threshold;

    for (round = 0; round < ransac_rounds; round++) {
        /* Pick 5 random points */
		dvec2 r_pts_inner[5], l_pts_inner[5];
        int indices[5];
        int num_hyp;
        double E[90];
        int inliers_hyp[10];
        int first_hyp = -1, first_hyp_idx = -1, second_hyp = -1;
        int best = 0;
        int num_ident = 0;
        int inliers = 0;

        choose(n, 5, indices);

        for (i = 0; i < 5; i++) {
            r_pts_inner[i] = r_pts_norm[indices[i]];
            l_pts_inner[i] = l_pts_norm[indices[i]];

            /* Check for degeneracy */
            if (r_pts_inner[i].x == l_pts_inner[i].x &&
                r_pts_inner[i].y == l_pts_inner[i].y)
                num_ident++;
        }
        
        if (num_ident >= 3)
            continue;  /* choose another 5 */
        
        generate_Ematrix_hypotheses(5, r_pts_inner, l_pts_inner, &num_hyp, E);
        
        for (i = 0; i < num_hyp; i++) {
            int best_inlier;
            double score = 0.0;

            double E2[9], tmp[9], F[9];
            memcpy(E2, E + 9 * i, 9 * sizeof(double));
            E2[0] = -E2[0];
            E2[1] = -E2[1];
            E2[3] = -E2[3];
            E2[4] = -E2[4];
            E2[8] = -E2[8];

            matrix_transpose_product(3, 3, 3, 3, K2_inv, E2, tmp);
            matrix_product(3, 3, 3, 3, tmp, K1_inv, F);

            inliers = evaluate_Ematrix(n, r_pts, l_pts, // r_pts_norm, l_pts_norm, 
                                       thresh_norm, F, // E + 9 * i, 
                                       &best_inlier, &score);
       
            if (inliers > max_inliers ||
                (inliers == max_inliers && score < min_score)) {
                best = 1;
                max_inliers = inliers;
                min_score = score;
                memcpy(E_best, E + 9 * i, sizeof(double) * 9);
                r_best = r_pts_norm[best_inlier];
                l_best = l_pts_norm[best_inlier];
            }

            inliers_hyp[i] = inliers;
        }

        if (best) {
            for (i = 0; i < num_hyp; i++) {
                if (inliers_hyp[i] > first_hyp) {
                    first_hyp = inliers_hyp[i];
                    first_hyp_idx = i;
                }
            }

            for (i = 0; i < num_hyp; i++) {
                if (i != first_hyp_idx && inliers_hyp[i] > second_hyp) {
                    second_hyp = inliers_hyp[i];
                }
            }

            // printf("first: %d, second: %d\n", first_hyp, second_hyp);
        }
    }

    if (max_inliers > 0) {
        int best_inlier;
        double score;
        // find_extrinsics_essential(E_best, r_best, l_best, R_out, t_out);
        int success = 
            find_extrinsics_essential_multipt(E_best, n, 
                                              r_pts_norm, l_pts_norm, 
                                              R_out, t_out);
        int inliers = 0;
        double E2[9];
        double tmp[9], F_best[9];

        if (success == 0) {
            free(r_pts_norm);
            free(l_pts_norm);
            return 0;
        }

        memcpy(E2, E_best, 9 * sizeof(double));
        E2[0] = -E2[0];
        E2[1] = -E2[1];
        E2[3] = -E2[3];
        E2[4] = -E2[4];
        E2[8] = -E2[8];

        matrix_transpose_product(3, 3, 3, 3, K2_inv, E2, tmp);
        matrix_product(3, 3, 3, 3, tmp, K1_inv, F_best);

        inliers = evaluate_Ematrix(n, r_pts, l_pts, // r_pts_norm, l_pts_norm, 
                                   thresh_norm, F_best, &best_inlier,
                                   &score);

        // matrix_print(3, 3, F_best);

        // printf("  inliers: %d / %d [score: %0.3e]\n", inliers, n, score);
        fflush(stdout);
    }

    // matrix_print(3, 3, E_best);

    delete[] r_pts_norm;
	delete[] l_pts_norm;

    return max_inliers;
}

