package gov.nist.math.jampack;

/**
 * Zmat implements general complex matrix stored in a rectangular array class Z.
 * 
 * @version Pre-alpha, 1999-02-24
 * @author G. W. Stewart
 */
public final class Zmat {

    boolean isPosSemiDefinite = false;

    /** The real part of the matrix */
    final double re[][];

    /** The imaginary part of the matrix */
    final double im[][];

    /** The number of rows */
    public final int nr;

    /** The number of columns */
    public final int nc;

    /**
     * Creates a Zmat and initializes its real and imaginary parts to a pair of
     * arrays.
     * 
     * @param re
     *            Contains the real part.
     * @param im
     *            Contains the imaginary part.
     * @exception ZException
     *                if the dimensions of re and im do not match
     */
    public Zmat(double re[][], double im[][]) throws ZException {
        nr = re.length;
        nc = re[0].length;
        if (nr != im.length || nc != im[0].length) {
            throw new ZException("Inconsistent array dimensions");
        }
        this.re = new double[nr][nc];
        this.im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                this.re[i][j] = re[i][j];
                this.im[i][j] = im[i][j];
            }
        }
    }

    /**
     * Creates a Zmat and initializes it to an array of class Z.
     */
    public Zmat(Z A[][]) {
        nr = A.length;
        nc = A[0].length;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                re[i][j] = A[i][j].re;
                im[i][j] = A[i][j].im;
            }
        }
    }

    /**
     * Creates a Zmat and initializes its real part to to an array of class
     * double. The imaginary part is set to zero.
     */
    public Zmat(double A[][]) {
        nr = A.length;
        nc = A[0].length;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                re[i][j] = A[i][j];
                im[i][j] = 0.0;
            }
        }
    }

    /**
     * Creates a Zmat and intitializes it to a Zmat.
     */
    public Zmat(Zmat A) {
        nr = A.nr;
        nc = A.nc;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                re[i][j] = A.re[i][j];
                im[i][j] = A.im[i][j];
            }
        }
    }

    /**
     * Creates a Zmat and initialize it to a Z1.
     */
    public Zmat(Z1 A) {
        nr = A.n;
        nc = 1;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            re[i][0] = A.re[i];
            im[i][0] = A.im[i];
        }
    }

    /**
     * Creates a Zmat and initialize it to a Zdiagmat.
     */
    public Zmat(Zdiagmat D) {
        nr = D.order;
        nc = D.order;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            re[i][i] = D.re(i);
            im[i][i] = D.im(i);
        }
    }

    /**
     * Creates a Zmat and initializes it to zero.
     */
    public Zmat(int nrow, int ncol) {
        this(nrow, ncol, false);
    }

    Zmat(int nrow, int ncol, boolean isPosSemiDefinite) {
        this.isPosSemiDefinite = isPosSemiDefinite;
        nr = nrow;
        nc = ncol;
        re = new double[nr][nc];
        im = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                re[i][j] = 0.0;
                im[i][j] = 0.0;
            }
        }
    }

    /**
     * Returns a copy of the real part of a Zmat.
     */
    public double[][] getRe() {
        double[][] A = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                A[i][j] = re[i][j];
            }
        }
        return A;
    }

    /**
     * Returns a copy of the imaginary part of a Zmat.
     */
    public double[][] getIm() {
        double[][] A = new double[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                A[i][j] = im[i][j];
            }
        }
        return A;
    }

    /**
     * Returns a copy of the real and imaginary parts as a complex array.
     */
    public Z[][] getZ() {
        Z[][] A = new Z[nr][nc];
        for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
                A[i][j] = new Z(re[i][j], im[i][j]);
            }
        }
        return A;
    }

    /**
     * Returns the (ii,jj)-element of a Zmat.
     * 
     * @param ii
     *            The row index of the element
     * @param jj
     *            The column index of the element
     */
    public Z get(int ii, int jj) {
        return new Z(re(ii - 1, jj - 1), im(ii - 1, jj - 1));
    }

    /**
     * Returns the zero-based (i,j)-element of a Zmat.
     * 
     * @param i
     *            The row index of the element
     * @param j
     *            The column index of the element
     */
    Z get0(int i, int j) {
        return new Z(re(i, j), im(i, j));
    }

    /**
     * Writes the (ii,jj) element of a Zmat.
     * 
     * @param ii
     *            The row index of the element
     * @param jj
     *            The column index of the element
     * @param a
     *            The new value of the element
     */
    public void put(int ii, int jj, Z a) {
        re[ii - 1][jj - 1] = a.re;
        im[ii - 1][jj - 1] = a.im;
    }

    /**
     * Writes the zero-based (i,j)-element of a Zmat.
     * 
     * @param i
     *            The row index of the element
     * @param j
     *            The column index of the element
     * @param a
     *            The new value of the element
     */
    void put0(int i, int j, Z a) {
        re[i][j] = a.re;
        im[i][j] = a.im;
    }

    /**
     * Returns the submatrix (ii1:ii2, jj1:jj2).
     * 
     * @param ii1
     *            The lower column index
     * @param ii2
     *            The upper column index
     * @param jj1
     *            The lower row index
     * @param jj2
     *            The upper row index
     */
    public Zmat get(int ii1, int ii2, int jj1, int jj2) {
        int nrow = ii2 - ii1 + 1;
        int ncol = jj2 - jj1 + 1;
        Zmat A = new Zmat(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                A.re[i][j] = re[i + ii1 - 1][j + jj1 - 1];
                A.im[i][j] = im[i + ii1 - 1][j + jj1 - 1];
            }
        }
        return A;
    }

    /**
     * Overwrites the submatrix (ii1:ii2, jj1:jj2) with a Zmat.
     * 
     * @param ii1
     *            The lower column index
     * @param ii2
     *            The upper column index
     * @param jj1
     *            The lower row index
     * @param jj2
     *            The upper row index
     * @param A
     *            The new value of the submatrix
     */
    public void put(int ii1, int ii2, int jj1, int jj2, Zmat A) {
        int nrow = ii2 - ii1 + 1;
        int ncol = jj2 - jj1 + 1;
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                re[i + ii1 - 1][j + jj1 - 1] = A.re[i][j];
                im[i + ii1 - 1][j + jj1 - 1] = A.im[i][j];
            }
        }
    }

    /**
     * Returns the submatrix (ii[], jj1:jj2).
     * 
     * @param ii
     *            Contains the row indices of the submatrix
     * @param jj1
     *            The lower column index
     * @param jj2
     *            The upper column index
     */
    public Zmat get(int ii[], int jj1, int jj2) {
        int nrow = ii.length;
        int ncol = jj2 - jj1 + 1;
        Zmat A = new Zmat(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                A.re[i][j] = re[ii[i] - 1][j + jj1 - 1];
                A.im[i][j] = im[ii[i] - 1][j + jj1 - 1];
            }
        }
        return A;
    }

    /**
     * Overwrites the submatrix (ii[], jj1:jj2) with a Zmat.
     * 
     * @param ii
     *            Contains the row indices of the submatrix
     * @param jj1
     *            The lower column index
     * @param jj2
     *            The upper column index
     * @param A
     *            The new value of the submatrix.
     */
    public void put(int ii[], int jj1, int jj2, Zmat A) {
        int nrow = ii.length;
        int ncol = jj2 - jj1 + 1;
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                re[ii[i] - 1][j + jj1 - 1] = A.re[i][j];
                im[ii[i] - 1][j + jj1 - 1] = A.im[i][j];
            }
        }
    }

    /**
     * Returns the submatrix (ii1:ii2, jj[]).
     * 
     * @param ii1
     *            The lower row index
     * @param ii2
     *            The upper row index
     * @param jj
     *            Contains the column indices of the submatrix
     */
    public Zmat get(int ii1, int ii2, int jj[]) {
        int nrow = ii2 - ii1 + 1;
        int ncol = jj.length;
        Zmat A = new Zmat(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                A.re[i][j] = re[i + ii1 - 1][jj[j] - 1];
                A.im[i][j] = im[i + ii1 - 1][jj[j] - 1];
            }
        }
        return A;
    }

    /**
     * Overwrites the submatrix (ii1:ii2, jj[]) with a Zmat.
     * 
     * @param ii1
     *            The lower row index
     * @param ii2
     *            The upper row index
     * @param jj
     *            Contains the column indices of the submatrix
     * @param A
     *            The new value of the submatrix
     */
    public void put(int ii1, int ii2, int jj[], Zmat A) {
        int nrow = ii2 - ii1 + 1;
        int ncol = jj.length;
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                re[i + ii1 - 1][jj[j] - 1] = A.re[i][j];
                im[i + ii1 - 1][jj[j] - 1] = A.im[i][j];
            }
        }
    }

    /**
     * Returns the submatrix (ii[], jj[]).
     * 
     * @param ii
     *            Contains the row indices of the submatrix
     * @param jj
     *            Contains the column indices of the submatrix
     */
    public Zmat get(int ii[], int jj[]) {
        int nrow = ii.length;
        int ncol = jj.length;
        Zmat A = new Zmat(nrow, ncol);
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                A.re[i][j] = re[ii[i] - 1][jj[j] - 1];
                A.im[i][j] = im[ii[i] - 1][jj[j] - 1];
            }
        }
        return A;
    }

    /**
     * Overwrites the submatrix (ii[], jj[]) with a Zmat. Returns the submatrix
     * (ii[], jj[])
     * 
     * @param ii
     *            Contains the row indices of the submatrix
     * @param jj
     *            Contains the column indices of the submatrix
     * @param A
     *            The value of the new submatrix
     */
    public void put(int ii[], int jj[], Zmat A) {
        int nrow = ii.length;
        int ncol = jj.length;
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
                re[ii[i] - 1][jj[j] - 1] = A.re[i][j];
                im[ii[i] - 1][jj[j] - 1] = A.im[i][j];
            }
        }
    }

    public double re(int i, int j) {
        return re[i][j];
    }

    public double im(int i, int j) {
        return im[i][j];
    }

    public void scaleRe(int i, int j, double scale) {
        re[i][j] *= scale;
    }

    public void scaleIm(int i, int j, double scale) {
        im[i][j] *= scale;
    }

    public void setRe(int i, int j, double val) {
        re[i][j] = val;
    }

    public void setIm(int i, int j, double val) {
        im[i][j] = val;
    }

    public int rows() {
        return nr;
    }

    public int cols() {
        return nc;
    }
}
