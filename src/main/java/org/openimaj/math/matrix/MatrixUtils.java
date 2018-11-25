/**
 * Copyright (c) 2011, The University of Southampton and the individual contributors.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *   *  Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *   *  Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 *   *  Neither the name of the University of Southampton nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package org.openimaj.math.matrix;

import gov.nist.math.jama.EigenvalueDecomposition;
import gov.nist.math.jama.JamaMatrix;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Random;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;
import ch.akuhn.matrix.SparseMatrix;

/**
 * Miscellaneous matrix operations.
 * 
 * @author Sina Samangooei (ss@ecs.soton.ac.uk)
 * @author Jonathon Hare (jsh2@ecs.soton.ac.uk)
 */
public final class MatrixUtils {

    private MatrixUtils() {
    }

    /**
     * Are any values NaN or Inf?
     * 
     * @param matrix
     *            matrix to test
     * @return true if any elements are NaN or Inf; false otherwise
     */
    public static boolean anyNaNorInf(JamaMatrix matrix) {
        for (double[] arrLine : matrix.getArray()) {
            for (double d : arrLine) {
                if (Double.isNaN(d) || Double.isInfinite(d)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Get the maximum absolute value of the diagonal.
     * 
     * @param matrix
     *            the matrix
     * @return the maximum absolute value
     */
    public static double maxAbsDiag(JamaMatrix matrix) {
        double max = -1.0;
        for (int i = 0; i < matrix.getColumnDimension(); i++) {
            double curr = Math.abs(matrix.get(i, i));
            if (max < curr) {
                max = curr;
            }
        }
        return max;
    }

    /**
     * Get the minimum absolute value of the diagonal.
     * 
     * @param matrix
     *            the matrix
     * @return the minimum absolute value
     */
    public static double minAbsDiag(JamaMatrix matrix) {
        double min = Double.MAX_VALUE;
        for (int i = 0; i < matrix.getColumnDimension(); i++) {
            double curr = Math.abs(matrix.get(i, i));
            if (min > curr) {
                min = curr;
            }
        }
        return min;
    }

    /**
     * Compute the principle square root, X, of the matrix A such that A=X*X
     * 
     * @param matrix
     *            the matrix
     * @return the sqrt of the matrix
     */
    public static JamaMatrix sqrt(JamaMatrix matrix) {
        // A = V*D*V'
        EigenvalueDecomposition evd = matrix.eig();
        JamaMatrix v = evd.getV();
        JamaMatrix d = evd.getD();

        // sqrt of cells of D and store in-place
        for (int r = 0; r < d.getRowDimension(); r++) {
            for (int c = 0; c < d.getColumnDimension(); c++) {
                d.set(r, c, Math.sqrt(d.get(r, c)));
            }
        }

        // Y = V*D/V
        // Y = V'.solve(V*D)'
        JamaMatrix a = v.inverse();
        JamaMatrix b = v.times(d).inverse();
        return a.solve(b).inverse();
    }

    /**
     * Computes the Moore-Penrose pseudoinverse. This is just a convenience
     * wrapper around {@link PseudoInverse#pseudoInverse(JamaMatrix)}.
     * 
     * @param matrix
     *            the matrix to invert.
     * @return the inverted matrix
     * @see PseudoInverse#pseudoInverse(JamaMatrix)
     */
    public static JamaMatrix pseudoInverse(JamaMatrix matrix) {
        return PseudoInverse.pseudoInverse(matrix);
    }

    /**
     * Compute the inverse square root, X, of the symmetric matrix A; A^-(1/2)
     * 
     * @param matrix
     *            the symmetric matrix
     * @return the inverse sqrt of the matrix
     */
    public static JamaMatrix invSqrtSym(JamaMatrix matrix) {
        // A = V*D*V'
        EigenvalueDecomposition evd = matrix.eig();
        JamaMatrix v = evd.getV();
        JamaMatrix d = evd.getD();

        // sqrt of cells of D and store in-place
        for (int r = 0; r < d.getRowDimension(); r++) {
            for (int c = 0; c < d.getColumnDimension(); c++) {
                if (d.get(r, c) > 0.0) {
                    d.set(r, c, 1.0 / Math.sqrt(d.get(r, c)));
                } else {
                    d.set(r, c, 0.0);
                }
            }
        }

        return v.times(d).times(v.transpose());
    }

    /**
     * Return a copy of the input matrix with all elements set to their absolute
     * value.
     * 
     * @param mat
     *            the matrix.
     * @return the absolute matrix.
     */
    public static JamaMatrix abs(JamaMatrix mat) {
        JamaMatrix copy = mat.copy();
        for (int i = 0; i < mat.getRowDimension(); i++) {
            for (int j = 0; j < mat.getColumnDimension(); j++) {
                copy.set(i, j, Math.abs(mat.get(i, j)));
            }
        }
        return copy;
    }

    /**
     * Check if two matrices are equal
     * 
     * @param m1
     *            first matrix
     * @param m2
     *            second matrix
     * @param eps
     *            epsilon for checking values
     * @return true if matrices have same size and all elements are equal within
     *         eps; false otherwise
     */
    public static boolean equals(JamaMatrix m1, JamaMatrix m2, double eps) {
        double[][] a1 = m1.getArray();
        double[][] a2 = m2.getArray();

        if (a1.length != a2.length || a1[0].length != a2[0].length) {
            return false;
        }

        for (int r = 0; r < a1.length; r++) {
            for (int c = 0; c < a1[r].length; c++) {
                if (Math.abs(a1[r][c] - a2[r][c]) > eps) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Return a copy of the matrix with all the values raised to a power.
     * 
     * @param mat
     *            the matrix.
     * @param exp
     *            the power.
     * @return a matrix.
     */
    public static JamaMatrix pow(JamaMatrix mat, double exp) {
        JamaMatrix copy = mat.copy();
        for (int i = 0; i < mat.getRowDimension(); i++) {
            for (int j = 0; j < mat.getColumnDimension(); j++) {
                copy.set(i, j, Math.pow(mat.get(i, j), exp));
            }
        }
        return copy;
    }

    /**
     * Generate a {@link String} representation of a matrix.
     * 
     * @param mat
     *            the matrix
     * @return a string representation
     */
    public static String toString(JamaMatrix mat) {
        StringWriter matWriter = new StringWriter();
        mat.print(new PrintWriter(matWriter), 5, 5);
        return matWriter.getBuffer().toString();
    }

    /**
     * Compute the sum of all elements of the matrix.
     * 
     * @param mat
     *            the matrix.
     * @return the sum.
     */
    public static double sum(JamaMatrix mat) {
        double sum = 0.0;
        for (int i = 0; i < mat.getRowDimension(); i++) {
            for (int j = 0; j < mat.getColumnDimension(); j++) {
                sum += mat.get(i, j);
            }
        }
        return sum;
    }

    /**
     * Zero the matrix.
     * 
     * @param m
     *            the matrix
     */
    public static void zero(JamaMatrix m) {
        m.timesEquals(0.0);
    }

    /**
     * Compute the real Eigen decomposition of a symmetric 2x2 matrix. Warning:
     * Doesn't check the size or whether the input is symmetric.
     * 
     * @param m
     *            the matrix
     * @return the Eigen vectors and values.
     */
    public static EigenValueVectorPair symmetricEig2x2(JamaMatrix m) {
        double a = m.get(0, 0);
        double b = m.get(0, 1);
        double c = b;
        double d = m.get(1, 1);

        double trace = a + d;
        double det = a * d - b * c;

        JamaMatrix val = new JamaMatrix(2, 2);
        double sqrtInner = (trace * trace / 4.0) - det;
        // FIXME: make it deal with imaginary numbers
        if (sqrtInner < 0.0) {
            EigenvalueDecomposition e = m.eig();
            return new EigenValueVectorPair(e.getD(), e.getV());
        }

        sqrtInner = Math.sqrt(sqrtInner);
        double firstEig = trace / 2.0 + sqrtInner;
        double secondEig = trace / 2.0 - sqrtInner;
        if (firstEig > secondEig) {
            double tmp = firstEig;
            firstEig = secondEig;
            secondEig = tmp;
        }

        val.set(0, 0, firstEig);
        val.set(1, 1, secondEig);

        JamaMatrix vec = new JamaMatrix(2, 2);

        double v1 = firstEig - a;
        double v2 = secondEig - a;
        double norm1 = Math.sqrt(v1 * v1 + b * b);
        double norm2 = Math.sqrt(b * b + v2 * v2);
        vec.set(0, 0, b / norm1);
        vec.set(0, 1, b / norm2);
        vec.set(1, 0, v1 / norm1);
        vec.set(1, 1, v2 / norm2);

        // to deal with rounding error
        vec.set(1, 0, vec.get(0, 1));

        return new EigenValueVectorPair(val, vec);
    }

    /**
     * An eigen decomposition that uses a deterministic method if the matrix is
     * 2x2.
     * 
     * This function returns values as in {@link EigenvalueDecomposition} i.e.
     * the largest eigen value is held in the [m.rows - 1,m.cols-1] (i.e. [1,1])
     * location
     * 
     * @param m
     * @return the decomposition
     */
    public static EigenValueVectorPair eig2x2(JamaMatrix m) {
        if (m.getColumnDimension() != 2 || m.getRowDimension() != 2) {
            EigenvalueDecomposition e = m.eig();
            return new EigenValueVectorPair(e.getD(), e.getV());
        }
        /**
         * A = 1 B = a + d C = ad - bc
         * 
         * x = ( - B (+/-) sqrt(B^2 - 4AC) )/ (2A)
         */
        double a = m.get(0, 0);
        double b = m.get(0, 1);
        double c = m.get(1, 0);
        double d = m.get(1, 1);

        double trace = a + d;
        double det = a * d - b * c;

        JamaMatrix val = new JamaMatrix(2, 2);
        double sqrtInner = (trace * trace / 4.0) - det;
        // FIXME: make it deal with imaginary numbers
        if (sqrtInner < 0.0) {
            EigenvalueDecomposition e = m.eig();
            return new EigenValueVectorPair(e.getD(), e.getV());
        }

        sqrtInner = Math.sqrt(sqrtInner);
        double firstEig = trace / 2.0 + sqrtInner;
        double secondEig = trace / 2.0 - sqrtInner;
        if (firstEig > secondEig) {
            double tmp = firstEig;
            firstEig = secondEig;
            secondEig = tmp;
        }

        val.set(0, 0, firstEig);
        val.set(1, 1, secondEig);

        JamaMatrix vec = new JamaMatrix(2, 2);
        if (b == 0.0 && c == 0.0) {
            vec.set(0, 0, 1.0);
            vec.set(1, 1, 1.0);
        } else {
            if (c != 0.0) {
                double v1 = firstEig - d;
                double v2 = secondEig - d;
                double norm1 = Math.sqrt(v1 * v1 + c * c);
                double norm2 = Math.sqrt(c * c + v2 * v2);
                vec.set(0, 0, v1 / norm1);
                vec.set(0, 1, v2 / norm2);
                vec.set(1, 0, c / norm1);
                vec.set(1, 1, c / norm2);
            } else if (b != 0.0) {
                double v1 = firstEig - a;
                double v2 = secondEig - a;
                double norm1 = Math.sqrt(v1 * v1 + b * b);
                double norm2 = Math.sqrt(b * b + v2 * v2);
                vec.set(0, 0, b / norm1);
                vec.set(0, 1, b / norm2);
                vec.set(1, 0, v1 / norm1);
                vec.set(1, 1, v2 / norm2);
            }
        }

        return new EigenValueVectorPair(val, vec);
    }

    /**
     * Construct a matrix from a 2D float array of data.
     * 
     * @param data
     *            the data.
     * @return the matrix.
     */
    public static JamaMatrix matrixFromFloat(float[][] data) {
        JamaMatrix out = new JamaMatrix(data.length, data[0].length);
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                out.set(j, i, data[i][j]);
            }
        }
        return out;
    }

    /**
     * Reduce the rank of a matrix by estimating the best (in a least-squares
     * sense) approximation using the thin SVD.
     * 
     * @param m
     *            the matrix to reduce.
     * @param rank
     *            the desired rank.
     * @return the rank-reduced matrix.
     */
    public static JamaMatrix reduceRank(JamaMatrix m, int rank) {
        if (rank > Math.min(m.getColumnDimension(), m.getRowDimension())) {
            return m;
        }

        DenseMatrix mjtA = new DenseMatrix(m.getArray());
        SVD svd;
        try {
            svd = SVD.factorize(mjtA);
        } catch (NotConvergedException e) {
            throw new RuntimeException(e);
        }

        DenseMatrix U = svd.getU();
        DenseMatrix Vt = svd.getVt();
        double[] svector = svd.getS();
        DenseMatrix S = new DenseMatrix(U.numColumns(), Vt.numRows());
        for (int i = 0; i < rank; i++) {
            S.set(i, i, svector[i]);
        }

        DenseMatrix C = new DenseMatrix(U.numRows(), S.numColumns());
        DenseMatrix out = new DenseMatrix(C.numRows(), Vt.numColumns());
        U.mult(S, C);
        C.mult(Vt, out);

        JamaMatrix outFinal = convert(out);
        return outFinal;
    }

    /**
     * Convert a MJT {@link DenseMatrix} to a {@link JamaMatrix}.
     * 
     * @param mjt
     *            {@link DenseMatrix} to convert
     * @return converted matrix.
     */
    public static JamaMatrix convert(DenseMatrix mjt) {
        return convert(mjt, mjt.numRows(), mjt.numColumns());
    }

    /**
     * Convert a MJT {@link DenseMatrix} to a {@link JamaMatrix}.
     * 
     * @param mjt
     *            {@link DenseMatrix} to convert
     * @param nrows
     *            number of rows to copy
     * @param ncols
     *            number of columns to copy
     * @return converted matrix.
     */
    public static JamaMatrix convert(DenseMatrix mjt, int nrows, int ncols) {
        double[][] d = new double[nrows][ncols];

        double[] mjtd = mjt.getData();
        for (int r = 0; r < nrows; r++) {
            for (int c = 0; c < ncols; c++) {
                d[r][c] = mjtd[r + c * d.length];
            }
        }
        return new JamaMatrix(d);
    }

    /**
     * Create a copy of a matrix with the columns in reverse order.
     * 
     * @param m
     *            the input matrix
     * @return a copy with the column order reversed
     */
    public static JamaMatrix reverseColumns(JamaMatrix m) {
        return reverseColumnsInplace(m.copy());
    }

    /**
     * Reverse the column order of the input matrix inplace.
     * 
     * @param m
     *            the input matrix
     * @return the input matrix
     */
    public static JamaMatrix reverseColumnsInplace(JamaMatrix m) {
        double[][] data = m.getArray();
        int rows = data.length;
        int cols = data[0].length;
        int halfCols = cols / 2;

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < halfCols; c++) {
                double tmp = data[r][c];
                data[r][c] = data[r][cols - c - 1];
                data[r][cols - c - 1] = tmp;
            }
        }

        return m;
    }

    /**
     * Create a copy of a matrix with the rows in reverse order.
     * 
     * @param m
     *            the input matrix
     * @return a copy with the row order reversed
     */
    public static JamaMatrix reverseRows(JamaMatrix m) {
        return reverseRowsInplace(m.copy());
    }

    /**
     * Reverse the row order of the input matrix inplace.
     * 
     * @param m
     *            the input matrix
     * @return the input matrix
     */
    public static JamaMatrix reverseRowsInplace(JamaMatrix m) {
        double[][] data = m.getArray();
        int rows = data.length;
        int halfRows = rows / 2;

        for (int r = 0; r < halfRows; r++) {
            double[] tmp = data[r];
            data[r] = data[rows - r - 1];
            data[rows - r - 1] = tmp;
        }

        return m;
    }

    /**
     * Create a diagonal matrix
     * 
     * @param s
     *            length diagonal numbers
     * @return new Matrix(s.length,s.length) s.t. diagonal element i,i = s[i]
     */
    public static JamaMatrix diag(double[] s) {
        JamaMatrix r = new JamaMatrix(s.length, s.length);
        for (int i = 0; i < s.length; i++) {
            r.set(i, i, s[i]);
        }
        return r;
    }

    /**
     * Set the values of the elements in a single column to a constant value.
     * 
     * @param m
     *            the matrix
     * @param c
     *            the column
     * @param v
     *            the constant value
     * @return the matrix
     */
    public static JamaMatrix setColumn(JamaMatrix m, int c, double v) {
        double[][] data = m.getArray();
        int rows = m.getRowDimension();

        for (int r = 0; r < rows; r++) {
            data[r][c] = v;
        }

        return m;
    }

    /**
     * Set the values of the elements in a single column to a constant value.
     * 
     * @param m
     *            the matrix
     * @param r
     *            the row
     * @param v
     *            the constant value
     * @return the matrix
     */
    public static JamaMatrix setRow(JamaMatrix m, int r, double v) {
        double[][] data = m.getArray();
        int cols = m.getColumnDimension();

        for (int c = 0; c < cols; c++) {
            data[r][c] = v;
        }

        return m;
    }

    /**
     * Fill a matrix with a constant value.
     * 
     * @param m
     *            the matrix
     * @param v
     *            the constant value
     * @return the matrix
     */
    public static JamaMatrix fill(JamaMatrix m, double v) {
        double[][] data = m.getArray();

        int rows = m.getRowDimension();
        int cols = m.getColumnDimension();

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                data[r][c] = v;
            }
        }

        return m;
    }

    /**
     * Subtract a constant from all values
     * 
     * @param m
     *            the matrix
     * @param v
     *            the constant value
     * @return the matrix
     */
    public static JamaMatrix minus(JamaMatrix m, double v) {
        double[][] data = m.getArray();

        int rows = m.getRowDimension();
        int cols = m.getColumnDimension();

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                data[r][c] -= v;
            }
        }

        return m;
    }

    /**
     * Add a constant to all values
     * 
     * @param m
     *            the matrix
     * @param v
     *            the constant value
     * @return the matrix
     */
    public static JamaMatrix plus(JamaMatrix m, double v) {
        double[][] data = m.getArray();

        int rows = m.getRowDimension();
        int cols = m.getColumnDimension();

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                data[r][c] += v;
            }
        }

        return m;
    }

    /**
     * Get a reshaped copy of the input matrix
     * 
     * @param m
     *            the matrix to reshape
     * @param newRows
     *            the new number of rows
     * @return new matrix
     */
    public static JamaMatrix reshape(JamaMatrix m, int newRows) {
        int oldCols = m.getColumnDimension();
        int length = oldCols * m.getRowDimension();
        int newCols = length / newRows;
        JamaMatrix mat = new JamaMatrix(newRows, newCols);

        double[][] m1v = m.getArray();
        double[][] m2v = mat.getArray();

        int r1 = 0, r2 = 0, c1 = 0, c2 = 0;
        for (int i = 0; i < length; i++) {
            m2v[r2][c2] = m1v[r1][c1];

            c1++;
            if (c1 >= oldCols) {
                c1 = 0;
                r1++;
            }

            c2++;
            if (c2 >= newCols) {
                c2 = 0;
                r2++;
            }
        }

        return mat;
    }

    /**
     * Get a reshaped copy of the input matrix
     * 
     * @param m
     *            the matrix to reshape
     * @param newRows
     *            the new number of rows
     * @param columnMajor
     *            if true, values are drawn and placed down columns first. If
     *            false values are drawn and placed across rows first
     * @return new matrix
     */
    public static JamaMatrix reshape(JamaMatrix m, int newRows, boolean columnMajor) {
        int oldCols = m.getColumnDimension();
        int oldRows = m.getRowDimension();
        int length = oldCols * m.getRowDimension();
        int newCols = length / newRows;
        JamaMatrix mat = new JamaMatrix(newRows, newCols);

        double[][] m1v = m.getArray();
        double[][] m2v = mat.getArray();

        int r1 = 0, r2 = 0, c1 = 0, c2 = 0;
        if (!columnMajor) {
            for (int i = 0; i < length; i++) {
                m2v[r2][c2] = m1v[r1][c1];

                c1++;
                if (c1 >= oldCols) {
                    c1 = 0;
                    r1++;
                }

                c2++;
                if (c2 >= newCols) {
                    c2 = 0;
                    r2++;
                }
            }
        }
        else {
            for (int i = 0; i < length; i++) {
                m2v[r2][c2] = m1v[r1][c1];

                r1++;
                if (r1 >= oldRows) {
                    r1 = 0;
                    c1++;
                }

                r2++;
                if (r2 >= newRows) {
                    r2 = 0;
                    c2++;
                }
            }
        }

        return mat;
    }

    /**
     * Compute the sum of values in a single column
     * 
     * @param m
     *            the matrix
     * @param col
     *            the column
     * @return the sum of values in column col
     */
    public static double sumColumn(JamaMatrix m, int col) {
        double[][] data = m.getArray();
        int rows = m.getRowDimension();

        double sum = 0.0;
        for (int r = 0; r < rows; r++) {
            sum += data[r][col];
        }

        return sum;
    }

    /**
     * Compute the sum of values in a single row
     * 
     * @param m
     *            the matrix
     * @param row
     *            the row
     * @return the sum of values in row row
     */
    public static double sumRow(JamaMatrix m, int row) {
        double[][] data = m.getArray();
        int cols = m.getColumnDimension();

        double sum = 0.0;
        for (int c = 0; c < cols; c++) {
            sum += data[row][c];
        }

        return sum;
    }

    /**
     * Increment values in a single column by a constant
     * 
     * @param m
     *            the matrix
     * @param col
     *            the column
     * @param value
     *            the constant
     * @return the matrix
     */
    public static JamaMatrix incrColumn(JamaMatrix m, int col, double value) {
        double[][] data = m.getArray();
        int rows = m.getRowDimension();

        for (int r = 0; r < rows; r++) {
            data[r][col] += value;
        }

        return m;
    }

    /**
     * Increment values in a single column by a constant
     * 
     * @param m
     *            the matrix
     * @param row
     *            the row
     * @param value
     *            the constant
     * @return the sum of values in row row
     */
    public static JamaMatrix incrRow(JamaMatrix m, int row, double value) {
        double[][] data = m.getArray();
        int cols = m.getColumnDimension();

        for (int c = 0; c < cols; c++) {
            data[row][c] += value;
        }

        return m;
    }

    /**
     * Round (using {@link Math#round(double)} each value of the matrix.
     * 
     * @param times
     * @return same matrix as handed in
     */
    public static JamaMatrix round(JamaMatrix times) {
        double[][] data = times.getArray();
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                data[i][j] = Math.round(data[i][j]);
            }
        }
        return times;
    }

    /**
     * min(A,B) returns an array the same size as A and B with the smallest
     * elements taken from A or B. The dimensions of A and B must match.
     * 
     * @param A
     * @param B
     * @return new Matrix filled with min from A and B
     */
    public static JamaMatrix min(JamaMatrix A, JamaMatrix B) {
        double[][] dataA = A.getArray();
        double[][] dataB = B.getArray();
        JamaMatrix ret = A.copy();
        double[][] dataRet = ret.getArray();
        for (int i = 0; i < dataA.length; i++) {
            for (int j = 0; j < dataB[i].length; j++) {
                dataRet[i][j] = Math.min(dataA[i][j], dataB[i][j]);
            }
        }
        return ret;
    }

    /**
     * d to the power of each value in range. As with {@link #range} range is: -
     * a single number (a) (0:1:a) - two numbers (a,b) (a:1:b) - three numbers
     * (a,b,c) (a:b:c)
     * 
     * Any other amount of range results in a {@link RuntimeException}.
     * 
     * @param d
     * @param range
     * @return d to the power of each value in range
     */
    public static JamaMatrix rangePow(double d, double... range) {
        double start, end, delta;
        if (range.length == 1) {
            start = 0.0;
            end = range[0];
            delta = 1.0;
        } else if (range.length == 2) {
            start = range[0];
            end = range[1];
            delta = 1.0;
        }
        else if (range.length == 3) {
            start = range[0];
            end = range[2];
            delta = range[1];
        }
        else {
            throw new RuntimeException("Invalid range options selected");
        }
        int l = (int) ((end - start + 1) / delta);
        double[][] out = new double[1][l];
        for (int i = 0; i < l; i++) {
            out[0][i] = Math.pow(d, start + (i * delta));
        }
        return new JamaMatrix(out);
    }

    /**
     * range is: - a single number (a) (0:1:a) - two numbers (a,b) (a:1:b) -
     * three numbers (a,b,c) (a:b:c)
     * 
     * @param range
     * @return the range defined
     */
    public static JamaMatrix range(double... range) {
        double start, end, delta;
        if (range.length == 1) {
            start = 0.0;
            end = range[0];
            delta = 1.0;
        } else if (range.length == 2) {
            start = range[0];
            end = range[1];
            delta = 1.0;
        }
        else if (range.length == 3) {
            start = range[0];
            end = range[2];
            delta = range[1];
        }
        else {
            throw new RuntimeException("Invalid range options selected");
        }
        int l = (int) Math.floor((end - start) / delta) + 1;
        double[][] out = new double[1][l];
        for (int i = 0; i < l; i++) {
            out[0][i] = start + (i * delta);
        }
        return new JamaMatrix(out);
    }

    /**
     * range is: - a single number (a) (0:1:a) - two numbers (a,b) (a:1:b) -
     * three numbers (a,b,c) (a:b:c)
     * 
     * @param range
     * @return the range defined
     */
    public static JamaMatrix range(int... range) {
        int start, end, delta;
        if (range.length == 1) {
            start = 0;
            end = range[0];
            delta = 1;
        } else if (range.length == 2) {
            start = range[0];
            end = range[1];
            delta = 1;
        }
        else if (range.length == 3) {
            start = range[0];
            end = range[2];
            delta = range[1];
        }
        else {
            throw new RuntimeException("Invalid range options selected");
        }
        int l = (int) Math.floor((end - start) / delta) + 1;
        double[][] out = new double[1][l];
        for (int i = 0; i < l; i++) {
            out[0][i] = start + (i * delta);
        }
        return new JamaMatrix(out);
    }

    /**
     * Given two row vectors, construct the power set of rowvector combinations.
     * 
     * @param A
     * @param B
     * @return a new matrix of size A.cols * B.cols
     */
    public static JamaMatrix ntuples(JamaMatrix A, JamaMatrix B) {
        double[][] Adata = A.getArray();
        double[][] Bdata = B.getArray();

        double[][] out = new double[2][Adata[0].length * Bdata[0].length];
        int i = 0;
        for (double a : Adata[0]) {
            for (double b : Bdata[0]) {
                out[0][i] = a;
                out[1][i] = b;
                i++;
            }
        }
        return new JamaMatrix(out);
    }

    /**
     * Given a matrix, repeat the matrix over i rows and j columns.
     * 
     * @param x
     * @param i
     * @param j
     * @return repeated matrix
     */
    public static JamaMatrix repmat(JamaMatrix x, int i, int j) {
        double[][] xdata = x.getArray();
        double[][] newmat = new double[xdata.length * i][xdata[0].length * j];
        for (int k = 0; k < newmat.length; k += xdata.length) {
            for (int l = 0; l < newmat[0].length; l += xdata[0].length) {
                int rowcopyindex = 0;
                for (double[] ds : xdata) {
                    System.arraycopy(ds, 0, newmat[k + rowcopyindex], l, xdata[0].length);
                    rowcopyindex += 1;
                }
            }
        }
        return new JamaMatrix(newmat);
    }

    /**
     * Horizontally stack all the matrices provided, i.e., ret = [x1 x2 x3 x4 ...
     * xn]
     * 
     * @param x
     * @return horizontally stacked
     */
    public static JamaMatrix hstack(JamaMatrix... x) {
        int height = x[0].getRowDimension();
        int width = 0;
        for (JamaMatrix matrix : x) {
            width += matrix.getColumnDimension();
        }
        double[][] newmat = new double[height][width];
        int colindex = 0;
        for (JamaMatrix matrix : x) {
            double[][] matdata = matrix.getArray();
            int w = matrix.getColumnDimension();
            for (int i = 0; i < height; i++) {
                System.arraycopy(matdata[i], 0, newmat[i], colindex, w);
            }
            colindex += w;
        }
        return new JamaMatrix(newmat);
    }

    /**
     * Add the rows to the mat at rowIndex. Assumes MANY things with no checks:
     * rows.rows == rowIndex.length, mat.cols == rows.cols, rowIndex.length &lt;
     * mat.rows, for x in rowIndex: x &lt; mat.rows &amp;&amp; x &gt;= 0 etc.
     * 
     * @param mat
     * @param rows
     * @param rowIndex
     * @return the input matrix
     */
    public static JamaMatrix plusEqualsRow(JamaMatrix mat, JamaMatrix rows, int[] rowIndex) {
        double[][] matdata = mat.getArray();
        double[][] rowdata = rows.getArray();
        int i = 0;
        for (int row : rowIndex) {
            for (int j = 0; j < rowdata[i].length; j++) {
                matdata[row][j] += rowdata[i][j];
            }
            i++;
        }
        return mat;
    }

    /**
     * @param x
     * @param val
     * @return a new matrix for x &lt; val
     */
    public static JamaMatrix lessThan(JamaMatrix x, double val) {
        JamaMatrix retMat = x.copy();
        double[][] data = x.getArray();
        double[][] retdata = retMat.getArray();
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                retdata[i][j] = data[i][j] < val ? 1.0 : 0.0;
            }
        }
        return retMat;
    }

    /**
     * @param x
     * @param val
     * @return a new matrix for x &gt; val
     */
    public static JamaMatrix greaterThan(JamaMatrix x, double val) {
        JamaMatrix retMat = x.copy();
        double[][] data = x.getArray();
        double[][] retdata = retMat.getArray();
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                retdata[i][j] = data[i][j] > val ? 1.0 : 0.0;
            }
        }
        return retMat;
    }

    /**
     * @param x
     * @return a new matrix for x1 &amp;&amp; x2 &amp;&amp; ... &amp;&amp; xn where &amp;&amp; means "!=0"
     */
    public static JamaMatrix and(JamaMatrix... x) {
        JamaMatrix retMat = MatrixUtils.ones(x[0].getRowDimension(), x[0].getColumnDimension());
        double[][] retdata = retMat.getArray();

        for (JamaMatrix matrix : x) {
            double[][] data = matrix.getArray();
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data[i].length; j++) {
                    retdata[i][j] = (retdata[i][j] != 0 && data[i][j] != 0) ? 1.0 : 0.0;
                }
            }
        }
        return retMat;
    }

    /**
     * @param rowDimension
     * @param columnDimension
     * @return matrix of dimensions filled with ones
     */
    public static JamaMatrix ones(int rowDimension, int columnDimension) {
        JamaMatrix ret = new JamaMatrix(rowDimension, columnDimension);
        return plus(ret, 1.0);
    }

    /**
     * @param x
     * @return logical-and each column of x
     */
    public static JamaMatrix all(JamaMatrix x) {
        int cols = x.getColumnDimension();
        int rows = x.getRowDimension();
        JamaMatrix ret = new JamaMatrix(1, cols);
        double[][] retdata = ret.getArray();
        double[][] data = x.getArray();
        for (int i = 0; i < cols; i++) {
            boolean cool = true;
            for (int j = 0; j < rows; j++) {
                cool = data[j][i] != 0 && cool;
                if (!cool) {
                    break;
                }
            }
            retdata[0][i] = cool ? 1.0 : 0.0;
        }
        return ret;
    }

    /**
     * @param vals
     * @return given vals, return the array indexes where vals != 0
     */
    public static int[] valsToIndex(double[] vals) {
        int nindex = 0;
        for (double d : vals) {
            nindex += (d != 0.0) ? 1 : 0;
        }

        int[] indexes = new int[nindex];
        nindex = 0;
        int i = 0;
        for (double d : vals) {
            if (d != 0.0) {
                indexes[i] = nindex;
                i++;
            }
            nindex++;
        }
        return indexes;
    }

    /**
     * For every value in x greater than {@code val} set {@code newVal}.
     * 
     * @param x
     * @param val
     * @param newVal
     * @return same matrix handed in
     */
    public static JamaMatrix greaterThanSet(JamaMatrix x, int val, int newVal) {
        double[][] data = x.getArray();
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                data[i][j] = data[i][j] > val ? newVal : data[i][j];
            }
        }
        return x;
    }

    /**
     * For every value in x less than {@code val} set {@code newVal}.
     * 
     * @param x
     * @param val
     * @param newVal
     * @return same matrix handed in
     */
    public static JamaMatrix lessThanSet(JamaMatrix x, int val, int newVal) {
        double[][] data = x.getArray();
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                data[i][j] = data[i][j] < val ? newVal : data[i][j];
            }
        }
        return x;
    }

    /**
     * Subtract the given row vector from every row of the given matrix,
     * returning the result in a new matrix.
     * 
     * @param in
     *            the matrix
     * @param row
     *            the row vector
     * @return the resultant matrix
     */
    public static JamaMatrix minusRow(JamaMatrix in, double[] row) {
        JamaMatrix out = in.copy();
        double[][] outData = out.getArray();
        int rows = out.getRowDimension();
        int cols = out.getColumnDimension();

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                outData[r][c] -= row[c];
            }
        }

        return out;
    }

    /**
     * Subtract the given col vector (held as a Matrix) from every col of the
     * given matrix, returning the result in a new matrix.
     * 
     * @param in
     *            the matrix
     * @param col
     *            the col Matrix (Only the first column is used)
     * @return the resultant matrix
     */
    public static JamaMatrix minusCol(JamaMatrix in, JamaMatrix col) {
        JamaMatrix out = in.copy();
        double[][] outData = out.getArray();
        int rows = out.getRowDimension();
        int cols = out.getColumnDimension();
        double[][] colArr = col.getArray();

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                outData[r][c] -= colArr[r][0];
            }
        }

        return out;
    }

    /**
     * Add a matrix to another inline.
     * 
     * @param result
     *            the matrix to add to
     * @param add
     *            the matrix to add
     * @return the result matrix
     */
    public static JamaMatrix plusEquals(JamaMatrix result, JamaMatrix add) {
        int rows = result.getRowDimension();
        int cols = result.getColumnDimension();

        double[][] resultData = result.getArray();
        double[][] addData = add.getArray();

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                resultData[r][c] += addData[r][c];
            }
        }

        return result;
    }

    /**
     * Multiply a matrix by a constant inplace, returning the matrix.
     * 
     * @param m
     *            the matrix
     * @param val
     *            the value to multiply by
     * @return the matrix
     */
    public static JamaMatrix times(JamaMatrix m, double val) {
        double[][] data = m.getArray();

        int rows = m.getRowDimension();
        int cols = m.getColumnDimension();

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                data[r][c] *= val;
            }
        }

        return m;
    }

    /**
     * Convert an MTJ matrix into a 2d double array.
     * 
     * @param mat
     * @return a double array
     */
    public static double[][] mtjToDoubleArray(DenseMatrix mat) {
        double[][] out = new double[mat.numRows()][mat.numColumns()];
        double[] data = mat.getData();
        for (int r = 0; r < out.length; r++) {
            double[] outr = out[r];
            for (int c = 0; c < out[0].length; c++) {
                outr[c] = data[r + c * out.length];
            }
        }
        return out;
    }

    /**
     * Compute the sum of values in all rows.
     * 
     * @param m
     *            the matrix
     * @return the sum of values across all cols in all rows
     */
    public static JamaMatrix sumRows(JamaMatrix m) {
        double[][] data = m.getArray();
        int cols = m.getColumnDimension();
        int rows = m.getRowDimension();

        JamaMatrix sum = new JamaMatrix(rows, 1);
        double[][] sumArr = sum.getArray();
        for (int c = 0; c < cols; c++) {
            for (int r = 0; r < rows; r++)
            {
                sumArr[r][0] += data[r][c];
            }
        }

        return sum;
    }

    /**
     * Compute the sum of values in all cols.
     * 
     * @param m
     *            the matrix
     * @return the sum of values across all rows in all cols
     */
    public static JamaMatrix sumCols(JamaMatrix m) {
        double[][] data = m.getArray();
        int cols = m.getColumnDimension();
        int rows = m.getRowDimension();

        JamaMatrix sum = new JamaMatrix(1, cols);
        double[][] sumArr = sum.getArray();
        for (int c = 0; c < cols; c++) {
            for (int r = 0; r < rows; r++)
            {
                sumArr[0][c] += data[r][c];
            }
        }

        return sum;
    }

    /**
     * Generate a matrix with Gaussian distributed randoms.
     * 
     * @param rows
     *            the number of rows
     * @param cols
     *            the number of columns
     * @return a matrix containing values drawn from a 0 mean 1.0 sdev gaussian
     */
    public static JamaMatrix randGaussian(int rows, int cols) {
        JamaMatrix m = new JamaMatrix(rows, cols);
        double[][] d = m.getArray();
        for (int row = 0; row < d.length; row++) {
            for (int col = 0; col < d[row].length; col++) {
                d[row][col] = r.nextGaussian();
            }
        }
        return m;
    }

    /**
     * Compute the sparsity (i.e. ratio of non-zero elements to matrix size) of
     * the given matrix.
     * 
     * @param matrix
     *            the matrix
     * @return the sparsity
     */
    public static double sparsity(SparseMatrix matrix) {
        double density = matrix.used() / ((double) matrix.rowCount() * (double) matrix.columnCount());
        return 1.0 - density;
    }

    /**
     * Extract the diagonal component from the given matrix.
     * 
     * @param cv
     *            the matrix
     * @return a new matrix with the diagonal values from the input matrix and
     *         other values set to zero.
     */
    public static JamaMatrix diag(JamaMatrix cv) {
        JamaMatrix d = new JamaMatrix(cv.getRowDimension(), cv.getColumnDimension());

        for (int i = 0; i < Math.min(cv.getRowDimension(), cv.getColumnDimension()); i++) {
            d.set(i, i, cv.get(i, i));
        }

        return d;
    }

    /**
     * Extract the diagonal component from the given matrix.
     * 
     * @param cv
     *            the matrix
     * @return a new matrix with the diagonal values from the input matrix and
     *         other values set to zero.
     */
    public static double[] diagVector(JamaMatrix cv) {
        double[] d = new double[Math.min(cv.getRowDimension(), cv.getColumnDimension())];

        for (int i = 0; i < Math.min(cv.getRowDimension(), cv.getColumnDimension()); i++) {
            d[i] = cv.get(i, i);
        }

        return d;
    }

    /**
     * Format a matrix as a single-line string suitable for using in MATLAB or
     * Octave.
     * 
     * @param mat
     *            the matrix to format
     * @return the string
     */
    public static String toMatlabString(JamaMatrix mat) {
        return "[" + toString(mat).trim().replace("\n ", ";").replace(" ", ",") + "]";
    }

    /**
     * Format a matrix as a single-line string suitable for using in Python.
     * 
     * @param mat
     *            the matrix to format
     * @return the string
     */
    public static String toPythonString(JamaMatrix mat) {
        return "[[" + toString(mat).trim().replace("\n ", "][").replace(" ", ",") + "]]";
    }

    /**
     * @param mat
     * @return trace of the matrix
     */
    public static double trace(JamaMatrix mat) {
        double sum = 0.0;
        for (int i = 0; i < mat.getRowDimension(); i++) {
            sum += mat.get(i, i);
        }
        return sum;
    }

    /**
     * Solves the system <code>Ax = 0</code>, returning the vector x as an
     * array. Internally computes the least-squares solution using the SVD of
     * <code>A</code>.
     * 
     * @param A
     *            the matrix describing the system
     * @return the solution vector
     */
    public static double[] solveHomogeneousSystem(JamaMatrix A) {
        return solveHomogeneousSystem(new DenseMatrix(A.getArray()));
    }

    /**
     * Solves the system <code>Ax = 0</code>, returning the vector x as an
     * array. Internally computes the least-squares solution using the SVD of
     * <code>A</code>.
     * 
     * @param A
     *            the matrix describing the system
     * @return the solution vector
     */
    public static double[] solveHomogeneousSystem(double[][] A) {
        return solveHomogeneousSystem(new DenseMatrix(A));
    }

    /**
     * Solves the system <code>Ax = 0</code>, returning the vector x as an
     * array. Internally computes the least-squares solution using the SVD of
     * <code>A</code>.
     * 
     * @param A
     *            the matrix describing the system
     * @return the solution vector
     */
    public static double[] solveHomogeneousSystem(DenseMatrix A) {
        try {
            SVD svd = SVD.factorize(A);

            double[] x = new double[svd.getVt().numRows()];
            int c = svd.getVt().numColumns() - 1;

            for (int i = 0; i < x.length; i++) {
                x[i] = svd.getVt().get(c, i);
            }

            return x;
        } catch (NotConvergedException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Format a matrix as a single-line string suitable for using in Java.
     * 
     * @param mat
     *            the matrix to format
     * @return the string
     */
    public static String toJavaString(JamaMatrix mat) {
        return "{" + toString(mat).trim().replaceAll("^", "{").replace("\n ", "},{").replace(" ", ",") + "}}";
    }

    /**
     * Construct a matrix from a row-packed (i.e. row-by-row) vector of the
     * data.
     * 
     * @param vector
     *            the row-packed vector
     * @param ncols
     *            the number of columns
     * @return the reconstructed matrix
     */
    public static JamaMatrix fromRowPacked(double[] vector, int ncols) {
        int nrows = vector.length / ncols;
        JamaMatrix a = new JamaMatrix(nrows, ncols);
        double[][] ad = a.getArray();
        for (int r = 0, i = 0; r < nrows; r++) {
            for (int c = 0; c < ncols; c++, i++) {
                ad[r][c] = vector[i];
            }
        }

        return a;
    }

    /**
     * Compute the covariance matrix of the given samples (assumed each sample
     * is a row).
     * 
     * @param m
     *            the samples matrix
     * @return the covariance matrix
     */
    public static JamaMatrix covariance(JamaMatrix m) {
        int N = m.getRowDimension();
        return times(m.transpose().times(m), 1.0 / (N > 1 ? N - 1 : N));
    }

    /**
     * For each element of X, sign(X) returns 1 if the element is greater than
     * zero, 0 if it equals zero and -1 if it is less than zero.
     * 
     * @param m
     *            the matrix
     * @return the sign matrix
     */
    public static JamaMatrix sign(JamaMatrix m) {
        JamaMatrix o = new JamaMatrix(m.getRowDimension(), m.getColumnDimension());
        double[][] md = m.getArray();
        double[][] od = o.getArray();

        for (int r = 0; r < o.getRowDimension(); r++) {
            for (int c = 0; c < o.getColumnDimension(); c++) {
                if (md[r][c] > 0) {
                    od[r][c] = 1;
                }
                if (md[r][c] < 0) {
                    od[r][c] = -1;
                }
            }
        }

        return o;
    }

    /**
     * Return a copy of the input matrix where every value is the exponential of
     * the elements, e to the X.
     * 
     * @param m
     *            the input matrix
     * @return the exponential matrix
     */
    public static JamaMatrix exp(JamaMatrix m) {
        JamaMatrix o = new JamaMatrix(m.getRowDimension(), m.getColumnDimension());
        double[][] md = m.getArray();
        double[][] od = o.getArray();

        for (int r = 0; r < o.getRowDimension(); r++) {
            for (int c = 0; c < o.getColumnDimension(); c++) {
                od[r][c] = Math.exp(md[r][c]);
            }
        }

        return o;
    }

    /**
     * Return a copy of the input matrix where every value is the hyperbolic
     * tangent of the elements.
     * 
     * @param m
     *            the input matrix
     * @return the tanh matrix
     */
    public static JamaMatrix tanh(JamaMatrix m) {
        JamaMatrix o = new JamaMatrix(m.getRowDimension(), m.getColumnDimension());
        double[][] md = m.getArray();
        double[][] od = o.getArray();

        for (int r = 0; r < o.getRowDimension(); r++) {
            for (int c = 0; c < o.getColumnDimension(); c++) {
                od[r][c] = Math.tanh(md[r][c]);
            }
        }

        return o;
    }

    private static final Random r = new Random();
    static {
        System.setProperty("com.github.fommil.netlib.LAPACK", "com.github.fommil.netlib.NativeRefLAPACK");
    }
}
