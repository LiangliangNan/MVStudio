#include <Eigen/Dense>
#include <cstdio>

extern "C" {

    void dgelsy_driver(double *A, double *b, double *x, int m, int n, int nrhs) {
        using namespace Eigen;

        if (m < n) {
            printf("Error: driver now only works when m >= n\n");
            return;
        }

        // Map row-major inputs
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> matA(A, m, n);
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> matB(b, m, nrhs);

        // 🔴 核心：带列主元 QR（等价 dgelsy）
        ColPivHouseholderQR<MatrixXd> qr(matA);

        // 求解
        MatrixXd matX = qr.solve(matB);   // n x nrhs

        // ⚠️ 可选：rank 检查（对应 LAPACK rank）
        int rank = qr.rank();
        // printf("Rank = %d\n", rank);

        // Copy to output x (row-major n x nrhs)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < nrhs; j++)
                x[i * nrhs + j] = matX(i, j);
    }

    void matrix_invert(int n, double *A, double *Ainv) {
        using namespace Eigen;

        assert(A != nullptr);
        assert(Ainv != nullptr);

        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> matA(A, n, n);

        // LU factorization (partial pivoting)
        PartialPivLU<MatrixXd> lu(matA);

        // 直接求逆（不要做 isInvertible 检查）
        MatrixXd matInv = lu.inverse();

        // Copy to output
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                Ainv[i * n + j] = matInv(i, j);
    }

    int dgeev_driver(int n, double *A, double *evec, double *eval) {
        using namespace Eigen;

        assert(A != nullptr);
        assert(evec != nullptr);
        assert(eval != nullptr);

        // Map row-major input
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> matA(A, n, n);

        // NaN check（保留原逻辑）
        for (int i = 0; i < n * n; i++) {
            if (A[i] != A[i]) {
                printf("[dgeev_driver] Error: nan encountered\n");
                return 0;
            }
        }

        // 🔴 EigenSolver（等价 dgeev）
        EigenSolver<MatrixXd> es(matA, /* computeEigenvectors = */ true);

        // ⚠️ Eigen 有 info()
        if (es.info() != Success) {
            printf("[dgeev_driver] Eigen decomposition failed\n");
            return 0;
        }

        VectorXcd eigenvalues = es.eigenvalues();     // complex
        MatrixXcd eigenvectors = es.eigenvectors();   // complex

        int count = 0;

        for (int i = 0; i < n; i++) {
            std::complex<double> lambda = eigenvalues(i);

            // 🔴 等价 wi[i] == 0
            if (std::abs(lambda.imag()) < 1e-12) {
                eval[count] = lambda.real();

                // 取对应特征向量（实部）
                for (int j = 0; j < n; j++) {
                    evec[count * n + j] = eigenvectors(j, i).real();
                }

                count++;
            }
        }

        return count;
    }


    void dgerqf_driver(int m, int n, double *A, double *R, double *Q)
    {
        using namespace Eigen;

        assert(A != nullptr);
        assert(R != nullptr);
        assert(Q != nullptr);

        // Map 输入（row-major）
        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> matA(A, m, n);

        // =========================
        // Step 1: A^T 做 QR
        // =========================
        MatrixXd AT = matA.transpose();   // n x m

        HouseholderQR<MatrixXd> qr(AT);

        // Q_t (n x n)
        MatrixXd Qt = qr.householderQ();

        // R_t (n x m) —— 只取上三角
        MatrixXd Rt = qr.matrixQR().topRows(n)
                            .template triangularView<Upper>();

        // =========================
        // Step 2: 构造 R = Rt^T
        // =========================
        MatrixXd Rmat = Rt.transpose();   // m x n

        // 强制成“右上三角”（完全对齐你原代码）
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (j < i)
                    R[i * n + j] = 0.0;
                else
                    R[i * n + j] = Rmat(i, j);
            }
        }

        // =========================
        // Step 3: 构造 Q = Qt^T
        // =========================
        MatrixXd Qmat = Qt.transpose();   // n x n

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Q[i * n + j] = Qmat(i, j);
            }
        }

        // =========================
        // （可选）数值验证（建议你调试时打开）
        // =========================
        /*
        MatrixXd Rm = Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(R, m, n);
        MatrixXd Qm = Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(Q, n, n);
        double err = (Rm * Qm - matA).norm();
        printf("RQ error = %e\n", err);
        */
    }


    int dgesvd_driver(int m, int n, double *A, double *U, double *S, double *VT)
    {
        using namespace Eigen;

        Map<Matrix<double, Dynamic, Dynamic, RowMajor>> matA(A, m, n);

        JacobiSVD<MatrixXd> svd(matA, ComputeFullU | ComputeFullV);

        VectorXd sing_vals = svd.singularValues();
        MatrixXd Ue = svd.matrixU();  // m x m
        MatrixXd Ve = svd.matrixV();  // n x n

        int k = std::min(m, n);
        for (int i = 0; i < k; i++) {
            S[i] = sing_vals(i);
        }

        // ✅ 保持原接口语义（最安全版本）
        // U = Ue
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                U[i * m + j] = Ue(i, j);
            }
        }

        // VT = V^T
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                VT[i * n + j] = Ve(j, i);
            }
        }

        return 1;
    }

    void dgesv_driver(int n, double *A, double *b, double *x) {
        assert(A && b && x);

        // 将 row-major C 数组转换为 Eigen 矩阵
        Eigen::MatrixXd matA(n, n);
        Eigen::VectorXd vecb(n);

        for (int i = 0; i < n; i++) {
            vecb(i) = b[i];
            for (int j = 0; j < n; j++) {
                matA(i, j) = A[i * n + j];
            }
        }

        // 求解线性系统 A x = b
        Eigen::FullPivLU<Eigen::MatrixXd> lu(matA);
        if (!lu.isInvertible()) {
            printf("[dgesv_driver] Warning: matrix is singular or nearly singular\n");
        }

        Eigen::VectorXd vecx = lu.solve(vecb);

        // 结果写回 C 数组
        for (int i = 0; i < n; i++) {
            x[i] = vecx(i);
        }
    }
}