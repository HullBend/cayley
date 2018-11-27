package ch.akuhn.matrix.eigenvalues;

import java.util.Arrays;

import com.github.fommil.netlib.ARPACK;
import org.netlib.util.doubleW;
import org.netlib.util.intW;

import ch.akuhn.matrix.KuhnMatrix;
import ch.akuhn.matrix.Vector;

/**
 * Finds a few eigenvalues of a matrix.
 * <p>
 * This class use ARPACK to find a few eigenvalues (&lambda;) and corresponding
 * eigenvectors (<b>x</b>) for the standard eigenvalue problem:
 *
 * <pre>
 * <code>A<b>x</b> = &lambda;<b>x</b></code>
 * </pre>
 *
 * where <code>A</code> is an <code>n</code> &times; <code>n</code> real
 * symmetric matrix.
 * <p>
 * The only thing that must be supplied in order to use this class on your
 * problem is to change the array dimensions appropriately, to specify
 * <em>which</em> eigenvalues you want to compute and to supply a matrix-vector
 * product
 *
 * <pre>
 * <b>w</b> &larr; A<b>v</b>
 * </pre>
 *
 * in the {@link #callback(Vector)} method.
 * <p>
 * Please refer to the ARPACK guide for further information.
 * <p>
 * <b>Example:</b>
 *
 * <pre>
 * Matrix A = <i>&hellip;square matrix&hellip;</i>;
 * Eigenvalues eigen = Eigenvalues.of(A).largest(4);
 * eigen.run();
 * double[] l = eigen.values;
 * Vector[] x = eigen.vectors;
 * </pre>
 *
 * @author Adrian Kuhn (Java) based on <code>ddsimp.f</code> by Richard Lehoucq,
 *         Danny Sorensen, Chao Yang (Fortran)
 *
 * @see "http://www.caam.rice.edu/software/ARPACK/UG"
 */
public abstract class FewEigenvalues extends Eigenvalues {

    private enum Which {
        LA, SA, LM, SM, BE
    };

    private Which which;

    public static FewEigenvalues of(final KuhnMatrix matrix) {
        if (!matrix.isSquare()) {
            throw new IllegalArgumentException("matrix is not square");
        }
        return new FewEigenvalues(matrix.columnCount()) {
            @Override
            protected Vector callback(Vector vector) {
                return matrix.mult(vector);
            }
        };
    }

    public FewEigenvalues(int n) {
        super(n);
        this.greatest(20);
    }

    private FewEigenvalues which(Which which, int nev) {
        this.which = which;
        this.nev = nev < n ? nev : n - 1;
        return this;
    }

    /** Compute the largest algebraic eigenvalues. */
    @Override
    public FewEigenvalues largest(int nev0) {
        return which(Which.LA, nev0);
    }

    /** Compute the smallest algebraic eigenvalues. */
    public FewEigenvalues smallest(int nev0) {
        return which(Which.SA, nev0);
    }

    /** Compute the largest eigenvalues in magnitude. */
    public FewEigenvalues greatest(int nev0) {
        return which(Which.LM, nev0);
    }

    /** Compute the smallest eigenvalues in magnitude. */
    public FewEigenvalues lowest(int nev0) {
        return which(Which.SM, nev0);
    }

    /**
     * Compute eigenvalues from both end of the spectrum. When the
     * <code>nev</code> is odd, compute one more from the high end than from the
     * low end.
     */
    public FewEigenvalues fromBothEnds(int nev0) {
        return which(Which.BE, nev0);
    }

    /**
     * Runs the eigenvalue decomposition, using an implicitly restarted Arnoldi
     * process (IRAP). Please refer to the ARPACK guide for more information.
     */
    @Override
    public Eigenvalues run() {
        if (n <= 0) {
            throw new IllegalArgumentException("n <= 0 : " + n);
        }
        if (nev <= 0) {
            throw new IllegalArgumentException("nev <= 0 : " + nev);
        }
        final ARPACK arpack = ARPACK.getInstance();
        /*
         * Setting up parameters for DSAUPD call.
         */
        final intW ido = new intW(0); // must start zero
        final String bmat = "I"; // standard problem
        final doubleW tol = new doubleW(0); // uses machine precision
        final intW info = new intW(0); // request random starting vector
        final double[] resid = new double[n]; // allocate starting vector
        /*
         * NVC is the largest number of basis vectors that will be used in the
         * Implicitly Restarted Arnoldi Process. Work per major iteration is
         * proportional to N*NCV*NCV.
         */
        final int ncv = Math.min(nev * 4, n); // rule of thumb use twice nev
        if (!(ncv > nev && ncv <= n)) {
            throw new IllegalArgumentException(
                    "!(ncv > nev && ncv <= n) : ncv = " + ncv + ", nev = " + nev + ", n = " + n);
        }
        final double[] v = new double[n * ncv];
        final double[] workd = new double[3 * n];
        final double[] workl = new double[ncv * (ncv + 8)];
        // Parameters
        final int[] iparam = new int[11];
        {
            final int ishfts = 1; // use the exact shift strategy
            final int maxitr = 300; // max number of Arnoldi iterations
            final int mode = 1; // use "mode 1"
            iparam[1 - 1] = ishfts;
            iparam[3 - 1] = maxitr;
            iparam[7 - 1] = mode;
        }
        final int[] ipntr = new int[11];
        // debug setting
        org.netlib.arpack.arpack_debug.ndigit.val = -3;
        org.netlib.arpack.arpack_debug.logfil.val = 6;
        org.netlib.arpack.arpack_debug.msgets.val = 0;
        org.netlib.arpack.arpack_debug.msaitr.val = 0;
        org.netlib.arpack.arpack_debug.msapps.val = 0;
        org.netlib.arpack.arpack_debug.msaupd.val = 0;
        org.netlib.arpack.arpack_debug.msaup2.val = 0;
        org.netlib.arpack.arpack_debug.mseigt.val = 0;
        org.netlib.arpack.arpack_debug.mseupd.val = 0;
        /*
         * Main loop (reverse communication loop)
         *
         * Repeatedly calls the routine DSAUPD and takes actions indicated by
         * parameter IDO until either convergence is indicated or maxitr has
         * been exceeded.
         */
        while (true) {
            arpack.dsaupd(
                    ido, // reverse communication parameter
                    bmat, // "I" = standard problem
                    n, // problem size
                    which.name(), // which values are requested?
                    nev, // how many values?
                    tol, // 0 = use machine precision
                    resid,
                    ncv,
                    v,
                    n,
                    iparam,
                    ipntr,
                    workd,
                    workl,
                    workl.length,
                    info);
            if (!(ido.val == 1 || ido.val == -1)) {
                break;
            }
            /*
             * Perform matrix vector multiplication
             *
             * y <--- OP*x
             *
             * The user should supply his own matrix-vector multiplication
             * routine here that takes workd(ipntr(1)) as the input, and return
             * the result to workd(ipntr(2)).
             */
            final int x0 = ipntr[1 - 1] - 1; // Fortran is off-by-one compared
                                             // to Java!
            final int y0 = ipntr[2 - 1] - 1;
            final Vector x = Vector.copy(workd, x0, n);
            final Vector y = this.callback(x);
            if (y.size() != n) {
                throw new IllegalArgumentException("callback.size() != n : " + y.size());
            }
            y.storeOn(workd, y0);
        }
        /*
         * Either we have convergence or there is an error.
         */
        if (info.val != 0) {
            throw new RuntimeException("dsaupd ERRNO = " + info.val
                    + ", see http://www.caam.rice.edu/software/ARPACK/UG/node136.html");
        }
        /*
         * Post-Process using DSEUPD.
         *
         * Computed eigenvalues may be extracted.
         *
         * The routine DSEUPD now called to do this post processing (other modes
         * may require more complicated post processing than "mode 1").
         */
        final boolean rvec = true; // request vectors
        final boolean[] select = new boolean[ncv];
        final double[] d = new double[ncv * 2];
        final double sigma = 0; // not used in "mode 1"
        final intW ierr = new intW(0);
        final intW nevW = new intW(nev);
        arpack.dseupd(
                rvec,
                "All",
                select,
                d,
                v,
                n,
                sigma,
                bmat,
                n,
                which.name(),
                nevW,
                tol.val,
                resid,
                ncv,
                v,
                n,
                iparam,
                ipntr,
                workd,
                workl,
                workl.length,
                ierr);
        if (ierr.val < 0) {
            throw new RuntimeException("dseupd ERRNO = " + info.val
                    + ", see http://www.caam.rice.edu/software/ARPACK/UG/node136.html");
        }
        /*
         * Eigenvalues are returned in the first column of the two dimensional
         * array D and the corresponding eigenvectors are returned in the first
         * NCONV (=IPARAM(5)) columns of the two dimensional array V if
         * requested. Otherwise, an orthogonal basis for the invariant subspace
         * corresponding to the eigenvalues in D is returned in V.
         */
        final int nconv = iparam[5 - 1];
        value = Arrays.copyOf(d, nconv);
        vector = new Vector[nconv];
        for (int i = 0; i < value.length; i++) {
            vector[i] = Vector.copy(v, i * n, n);
        }
        return this;
    }

    protected abstract Vector callback(Vector vector);
}
