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


}