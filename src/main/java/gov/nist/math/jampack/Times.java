package gov.nist.math.jampack;

/**
 * Times provides static methods to compute matrix products.
 * 
 * @version Pre-alpha, 1999-02-24
 * @author G. W. Stewart
 */
public final class Times {

    /**
     * Computes the product of a Z and a Zmat.
     * 
     * @param z
     *            The complex scalar
     * @param A
     *            The Zmat
     * @return zA
     */
    public static Zmat o(Z z, Zmat A) {
        Zmat B = new Zmat(A.nr, A.nc);
        for (int i = 0; i < A.nr; i++)
            for (int j = 0; j < A.nc; j++) {
                double reij = A.re(i, j);
                double imij = A.im(i, j);
                B.setRe(i, j, z.re * reij - z.im * imij);
                B.setIm(i, j, z.im * reij + z.re * imij);
            }
        return B;
    }

    /**
     * Computes the product of two Zmats.
     * 
     * @param A
     *            The first Zmat
     * @param B
     *            The second Zmat
     * @return AB
     * @exception ZException
     *                for unconformity
     */
    public static Zmat o(Zmat A, Zmat B) throws ZException {
        if (A.nc != B.nr) {
            throw new ZException("Unconformity in product");
        }
        Zmat C = new Zmat(A.nr, B.nc);
        for (int i = 0; i < A.nr; i++) {
            for (int k = 0; k < A.nc; k++) {
                for (int j = 0; j < B.nc; j++) {
                    double imkj = B.im(k, j);
                    double rekj = B.re(k, j);
                    double imik = A.im(i, k);
                    double reik = A.re(i, k);
                    C.setRe(i, j, C.re(i, j) + reik * rekj - imik * imkj);
                    C.setIm(i, j, C.im(i, j) + imik * rekj + reik * imkj);
                }
            }
        }
        return C;
    }

    /**
     * Computes A<sup>H</sup>A, where A is a Zmat.
     * 
     * @param A
     *            The Zmat
     * @return A<sup>H</sup>A
     */
    public static Zmat aha(Zmat A) {
        Zmat C = new Zmat(A.nc, A.nc, true);
        for (int k = 0; k < A.nr; k++) {
            for (int i = 0; i < A.nc; i++) {
                C.re[i][i] = C.re[i][i] + A.re[k][i] * A.re[k][i] + A.im[k][i] * A.im[k][i];
                C.im[i][i] = 0.;
                for (int j = i + 1; j < A.nc; j++) {
                    C.re[i][j] = C.re[i][j] + A.re[k][i] * A.re[k][j] + A.im[k][i] * A.im[k][j];
                    C.im[i][j] = C.im[i][j] + A.re[k][i] * A.im[k][j] - A.im[k][i] * A.re[k][j];
                }
            }
        }
        for (int i = 0; i < A.nc; i++) {
            for (int j = i + 1; j < A.nc; j++) {
                C.re[j][i] = C.re[i][j];
                C.im[j][i] = -C.im[i][j];
            }
        }
        return C;
    }

    /**
     * Computes AA<sup>H</sup>, where A is a Zmat.
     * 
     * @param A
     *            The Zmat
     * @return AA<sup>H</sup>
     */
    public static Zmat aah(Zmat A) {
        Zmat C = new Zmat(A.nr, A.nr, true);
        for (int i = 0; i < A.nr; i++) {
            for (int k = 0; k < A.nc; k++) {
                C.re[i][i] = C.re[i][i] + A.re[i][k] * A.re[i][k] + A.im[i][k] * A.im[i][k];
            }
            C.im[i][i] = 0.;
            for (int j = i + 1; j < A.nr; j++) {
                for (int k = 0; k < A.nc; k++) {
                    C.re[i][j] = C.re[i][j] + A.re[i][k] * A.re[j][k] + A.im[i][k] * A.im[j][k];
                    C.im[i][j] = C.im[i][j] - A.re[i][k] * A.im[j][k] + A.im[i][k] * A.re[j][k];
                }
                C.re[j][i] = C.re[i][j];
                C.im[j][i] = -C.im[i][j];
            }
        }
        return C;
    }

    /**
     * Computes the product of a Z and a Zdiagmat.
     * 
     * @param z
     *            The complex scalar
     * @param D
     *            The Zdiagmat
     * @return zD
     */
    public static Zdiagmat o(Z z, Zdiagmat D) {
        Zdiagmat B = new Zdiagmat(D);
        for (int i = 0; i < D.order; i++) {
            double rei = D.re(i);
            double imi = D.im(i);
            B.setRe(i, z.re * rei - z.im * imi);
            B.setIm(i, z.im * rei + z.re * imi);
        }
        return B;
    }

    /**
     * Computes the product of two Zdiagmats.
     * 
     * @param D1
     *            The first Zdiagmat
     * @param D2
     *            The second Zdiagmat
     * @return D1*D2
     * @exception ZException
     *                for unconformity
     */
    public static Zdiagmat o(Zdiagmat D1, Zdiagmat D2) throws ZException {
        if (D1.order != D2.order) {
            throw new ZException("Unconformity in product");
        }
        Zdiagmat D3 = new Zdiagmat(D1.order);
        for (int i = 0; i < D3.order; i++) {
            double d1rei = D1.re(i);
            double d1imi = D1.im(i);
            double d2rei = D2.re(i);
            double d2imi = D2.im(i);
            D3.setRe(i, d1rei * d2rei - d1imi * d2imi);
            D3.setIm(i, d1rei * d2imi + d1imi * d2rei);
        }
        return D3;
    }

    /**
     * Computes the product of a Zdiagmat and a Zmat.
     * 
     * @param D
     *            The Zdiagmat
     * @param A
     *            The Zmat
     * @return DA
     * @exception ZException
     *                for unconformity
     */
    public static Zmat o(Zdiagmat D, Zmat A) throws ZException {
        if (D.order != A.nr) {
            throw new ZException("Unconformity in product.");
        }
        Zmat B = new Zmat(A.nr, A.nc);
        for (int i = 0; i < A.nr; i++) {
            for (int j = 0; j < A.nc; j++) {
                double rei = D.re(i);
                double imi = D.im(i);
                double reij = A.re(i, j);
                double imij = A.im(i, j);
                B.setRe(i, j, rei * reij - imi * imij);
                B.setIm(i, j, rei * imij + imi * reij);
            }
        }
        return B;
    }

    /**
     * Computes the product of a Zmat and a Zdiagmat.
     * 
     * @param A
     *            The Zgmat
     * @param D
     *            The Zdiagmat
     * @return AD
     * @exception ZException
     *                for unconformity
     */
    public static Zmat o(Zmat A, Zdiagmat D) throws ZException {
        if (D.order != A.nc) {
            throw new ZException("Unconformity in product.");
        }
        Zmat B = new Zmat(A.nr, A.nc);
        for (int i = 0; i < A.nr; i++) {
            for (int j = 0; j < A.nc; j++) {
                double rej = D.re(j);
                double imj = D.im(j);
                double reij = A.re(i, j);
                double imij = A.im(i, j);
                B.setRe(i, j, rej * reij - imj * imij);
                B.setIm(i, j, rej * imij + imj * reij);
            }
        }
        return B;
    }
}
