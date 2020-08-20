package gov.nist.math.jampack;

/**
 * Zdiagmat is a storage efficient representation of a complex diagonal matrix.
 * 
 * @version Pre-alpha
 * @author G. W. Stewart
 */
public final class Zdiagmat {

    /** The order of the matrix */
    final int order;

    /** The real part of the diagonal */
    double re[];

    /** The imaginary part of the diagonal */
    final double im[];

    /** The order of the matrix (public) */
    public final int n;

    /** The index of the last diagonal (public) */
    final int dx;

    /**
     * Constructs a Zdiagmat and initializes it to zero.
     * 
     * @param order
     *            The order of the new Zdiagmat
     */
    public Zdiagmat(int order) {
        this.order = order;
        dx = order;
        n = order;
        re = new double[n];
        im = new double[n];
    }

    /**
     * Constructs a Zdiagmat and initializes it to a constant.
     * 
     * @param order
     *            The order of the new Zdiagmat
     * @param val
     *            The value to which the diagonal is to be initialized
     */
    public Zdiagmat(int order, Z val) {
        this.order = order;
        dx = order;
        n = order;
        re = new double[n];
        im = new double[n];
        for (int i = 0; i < n; i++) {
            re[i] = val.re;
            im[i] = val.im;
        }
    }

    /**
     * Constructs a Zdiagmat and initializes it to a Z1.
     * 
     * @param val a Z1
     */
    public Zdiagmat(Z1 val) {
        order = val.re.length;
        dx = order;
        n = order;
        re = new double[n];
        im = new double[n];
        for (int i = 0; i < n; i++) {
            re[i] = val.re[i];
            im[i] = val.im[i];
        }
    }

    /**
     * Constructs a Zdiagmat and initializes it to the diagonal of a Zmat.
     * 
     * @param A
     *            The Zmat
     * @param k
     *            The diagonal. For k=0 gives the princpal diagonal; k>0, the
     *            kth superdiagonal; k<0, the kth subdiagonal.
     * @exception ZException
     *                Thrown for k to large or small.
     */
    public Zdiagmat(Zmat A, int k) throws ZException {
        if (k >= 0) {
            if (k >= A.ncol) {
                throw new ZException("Diagonal out of range.");
            }
            order = Math.min(A.nrow, A.ncol - k);
            re = new double[order];
            im = new double[order];
            for (int i = 0; i < order; i++) {
                re[i] = A.re[i][i + k];
                im[i] = A.im[i][i + k];
            }
        } else {
            k = -k;
            if (k >= A.nrow) {
                throw new ZException("Diagonal out of range.");
            }
            order = Math.min(A.nrow - k, A.ncol);
            re = new double[order];
            im = new double[order];
            for (int i = 0; i < order; i++) {
                re[i] = A.re[i + k][i];
                im[i] = A.im[i + k][i];
            }
        }
        dx = order;
        n = order;
    }

    /**
     * Constructs a Zdiagmat and initializes it to the principal diagonal of a
     * Zmat.
     * 
     * @param A a Zmat
     * @exception ZException
     *                Passed from below.
     */
    public Zdiagmat(Zmat A) throws ZException {
        this(A, 0);
    }

    /**
     * Constructs a Zdiagmat and initializes it to another Zdiagmat.
     * 
     * @param D a Zdiagmat
     */
    public Zdiagmat(Zdiagmat D) {
        order = D.order;
        dx = order;
        n = order;
        re = new double[n];
        im = new double[n];

        for (int i = 0; i < n; i++) {
            re[i] = D.re[i];
            im[i] = D.im[i];
        }
    }

    /**
     * Gets the ii-th diagonal element of a Zdiagmat.
     * 
     * @param ii
     *            An integer
     * @return The ii-th element of this Zdiagmat
     */
    public Z get(int ii) {
        return new Z(re[ii - 1], im[ii - 1]);
    }

    /**
     * Writes the ii-th diagonal element of a Zdiagmat.
     * 
     * @param ii
     *            An integer
     * @param val a Z
     */
    public void put(int ii, Z val) {
        re[ii - 1] = val.re;
        im[ii - 1] = val.im;
    }

    public int order() {
        return order;
    }
}
