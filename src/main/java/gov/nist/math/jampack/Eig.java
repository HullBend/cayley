package gov.nist.math.jampack;

/**
 * Eig implements the eigenvalue-vector decomposition of of a square matrix.
 * Specifically given a diagonalizable matrix A, there is a non-singular matrix
 * X such that
 * 
 * <pre>
 *      D = X<sup>-1</sup> AX
 * </pre>
 * 
 * is diagonal. The columns of X are eigenvectors of A corresponding to the
 * diagonal elements of D. Eig implements X as a Zmat and D as a Zdiagmat.
 * <p>
 * Warning: if A is defective rounding error will allow Eig to compute a set of
 * eigenvectors. However, the matrix X will be ill conditioned.
 * 
 * @version Pre-alpha
 * @author G. W. Stewart
 */
public final class Eig {

    /** The matrix of eigenvectors */
    public Zmat X;

    /** The diagonal matrix of eigenvalues */
    public Zdiagmat D;

    /**
     * Creates an eigenvalue-vector decomposition of a square matrix A.
     * 
     * @param A
     *            The matrix whose decomposition is to be computed
     * @exception ZException
     *                Thrown if A is not square. <br>
     *                Passed from below.
     */
    public Eig(Zmat A) throws ZException {
        int i, j, k;
        double norm, scale;
        Z z, d;

        if (A.nr != A.nc) {
            throw new ZException("Matrix not square.");
        }

        int n = A.nr;

        /* Compute the Schur decomposition of $A$ and set up T and D. */
        Schur S = new Schur(A);

        Zmat T = S.T;

        D = new Zdiagmat(T);

        norm = Norm.fro(A);

        X = new Zmat(n, n);

        /* Compute the eigenvectors of T */
        for (k = n - 1; k >= 0; k--) {

            d = T.get0(k, k);

            X.re[k][k] = 1.0;
            X.im[k][k] = 0.0;

            for (i = k - 1; i >= 0; i--) {

                X.re[i][k] = -T.re[i][k];
                X.im[i][k] = -T.im[i][k];

                for (j = i + 1; j < k; j++) {

                    X.re[i][k] = X.re[i][k] - T.re[i][j] * X.re[j][k] + T.im[i][j] * X.im[j][k];
                    X.im[i][k] = X.im[i][k] - T.re[i][j] * X.im[j][k] - T.im[i][j] * X.re[j][k];
                }

                z = T.get0(i, i);
                z.minus(z, d);
                if (z.re == 0.0 && z.im == 0.0) { // perturb zero diagonal
                    z.re = 1.0e-16 * norm; // to avoid division by zero
                }
                z.div(X.get0(i, k), z);
                X.put0(i, k, z);
            }

            /* Scale the vector so its norm is one. */
            scale = 1.0 / Norm.fro(X, 1, X.rx, 1 + k, 1 + k);
            for (i = 0; i < X.nr; i++) {
                X.re[i][k] = scale * X.re[i][k];
                X.im[i][k] = scale * X.im[i][k];
            }
        }
        X = Times.o(S.U, X);
    }
}
