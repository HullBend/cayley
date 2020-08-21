package gov.nist.math.jampack;

/**
 * H computes the conjugate transpose of a matrix and the transpose of a complex
 * matrix.
 * 
 * @version Pre-alpha, 1999-02-24
 * @author G. W. Stewart
 */
public final class H {

    /**
     * Returns the conjugate transpose of a Zmat.
     * 
     * @param A
     *            The matrix to be conjugated and transposed
     * @return The conjugate transpose of A
     */
    public static Zmat o(Zmat A) {
        Zmat Ah = new Zmat(A.nc, A.nr);
        for (int i = 0; i < A.nr; i++)
            for (int j = 0; j < A.nc; j++) {
                Ah.setRe(j, i, A.re(i, j));
                Ah.setIm(j, i, -A.im(i, j));
            }
        return Ah;
    }

    /**
     * Returns the conjugate transpose of a Zdiagmat.
     * 
     * @param D
     *            The matrix to be conjugated and transposed
     * @return The conjugate transpose of D
     */
    public static Zdiagmat o(Zdiagmat D) {
        Zdiagmat Dh = new Zdiagmat(D);
        for (int i = 0; i < Dh.order; i++) {
            Dh.im[i] = -Dh.im[i];
        }
        return Dh;
    }

    /**
     * Returns the transpose of a Zmat.
     * 
     * @param A
     *            The matrix to be transposed
     * @return The transpose of A
     */
    public static Zmat trans(Zmat A) {
        Zmat Ah = new Zmat(A.nc, A.nr);
        for (int i = 0; i < A.nr; i++)
            for (int j = 0; j < A.nc; j++) {
                Ah.setRe(j, i, A.re(i, j));
                Ah.setIm(j, i, A.im(i, j));
            }
        return Ah;
    }
}
