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

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
import no.uib.cipr.matrix.sparse.SparseVector;

import ch.akuhn.matrix.KuhnDenseMatrix;
import ch.akuhn.matrix.DenseVector;
import ch.akuhn.matrix.KuhnMatrix;
import ch.akuhn.matrix.SparseMatrix;
import ch.akuhn.matrix.KuhnSparseVector;
import ch.akuhn.matrix.Vector;
import ch.akuhn.matrix.Vector.Entry;
import gov.nist.math.jama.JamaMatrix;

/**
 * Some helpful operations on {@link KuhnMatrix} instances from Adrian Kuhn's
 * library.
 *
 * @author Sina Samangooei (ss@ecs.soton.ac.uk)
 */
public final class MatlibMatrixUtils {

    private static final double EPS = 1e-8;

    /**
     * Compute the matrix sparsity (i.e. proportion of zero elements)
     *
     * @param mat
     *            the matrix
     * @return the sparsity
     * @see SparseMatrix#density()
     */
    public static double sparsity(KuhnMatrix mat) {
        return 1.0 - mat.density();
    }

    /**
     * Raise each element to the power d, operating on the matrix itself
     *
     * @param matrix
     *            the matrix
     * @param d
     *            the power
     * @return the input matrix, with elements raised to the given power
     */
    public static <T extends KuhnMatrix> T powInplace(T matrix, double d) {
        int rowN = 0;
        for (Vector ent : matrix.rows()) {
            for (Entry row : ent.entries()) {
                matrix.put(rowN, row.index, Math.pow(row.value, d));
            }
            rowN++;
        }
        return matrix;
    }

    /**
     * Left multiply two matrices: <code>R = D . A</code>
     *
     * @param D
     *            first matrix
     * @param A
     *            second matrix
     * @return result of multiplication
     */
    public static SparseMatrix times(DiagonalMatrix D, SparseMatrix A) {
        SparseMatrix mat = new SparseMatrix(A.rowCount(), A.columnCount());
        double[] Dvals = D.getVals();
        int rowIndex = 0;

        for (Vector row : A.rows()) {
            for (Entry ent : row.entries()) {
                mat.put(rowIndex, ent.index, ent.value * Dvals[rowIndex]);
            }
            rowIndex++;
        }

        return mat;
    }

    /**
     * Right multiply two matrices: <code>R = A . D</code>
     *
     * @param D
     *            first matrix
     * @param A
     *            second matrix
     * @return result of multiplication
     */
    public static SparseMatrix times(SparseMatrix A, DiagonalMatrix D) {
        SparseMatrix mat = new SparseMatrix(A.rowCount(), A.columnCount());
        int rowIndex = 0;
        double[] Dvals = D.getVals();
        for (Vector row : A.rows()) {
            for (Entry ent : row.entries()) {
                mat.put(rowIndex, ent.index, ent.value * Dvals[ent.index]);
            }
            rowIndex++;
        }
        return mat;
    }

    /**
     * Add two matrices, storing the results in the first:
     * <code>A = A + B</code>
     *
     * @param A
     *            first matrix
     * @param B
     *            matrix to add
     * @return A first matrix
     */
    public static SparseMatrix plusInplace(SparseMatrix A, KuhnMatrix B) {
        for (int i = 0; i < A.rowCount(); i++) {
            A.addToRow(i, B.row(i));
        }
        return A;
    }

    /**
     * Add two matrices, storing the results in the first:
     * <code>A = A + B</code>
     *
     * @param A
     *            first matrix
     * @param B
     *            matrix to add
     * @return A first matrix
     */
    @SuppressWarnings("unchecked")
    public static <T extends KuhnMatrix> T plusInplace(T A, KuhnMatrix B) {
        if (A instanceof SparseMatrix) {
            return (T) plusInplace((SparseMatrix) A, B);
        }
        for (int i = 0; i < A.rowCount(); i++) {
            Vector brow = B.row(i);
            for (int j = 0; j < A.columnCount(); j++) {
                A.row(i).add(j, brow.get(j));
            }
        }
        return A;
    }

    /**
     * Add a constant inplace <code>A = A + d</code>
     *
     * @param A
     *            first matrix
     * @param d
     *            the constant to add
     * @return A first matrix
     */
    public static <T extends KuhnMatrix> T plusInplace(T A, double d) {
        for (int i = 0; i < A.rowCount(); i++) {
            for (int j = 0; j < A.columnCount(); j++) {
                A.row(i).add(j, d);
            }
        }
        return A;
    }

    /**
     * Subtract two matrices, storing the result in the second:
     * <code>A = D - A</code>
     *
     * @param D
     *            first matrix
     * @param A
     *            second matrix
     * @return second matrix
     *
     */
    public static <T extends KuhnMatrix> T minusInplace(DiagonalMatrix D, T A) {
        double[] Dval = D.getVals();
        for (int i = 0; i < Dval.length; i++) {
            Iterable<Entry> rowents = A.row(i).entries();
            for (Entry entry : rowents) {
                A.put(i, entry.index, -entry.value);
            }
            A.put(i, i, Dval[i] - A.get(i, i));
        }
        return A;
    }

    /**
     * Subtract two matrices, storing the result in the first:
     * <code>A = A - B</code>
     *
     * @param A
     *            first matrix
     * @param B
     *            second matrix
     * @return first matrix
     *
     */
    public static KuhnMatrix minusInplace(KuhnMatrix A, KuhnMatrix B) {
        for (int i = 0; i < A.rowCount(); i++) {
            Iterable<Entry> rowents = A.row(i).entries();
            for (Entry entry : rowents) {
                A.put(i, entry.index, entry.value - B.get(i, entry.index));
            }
        }
        return A;
    }

    /**
     * Subtract a vector from another vector <code>A = A - D</code>
     *
     * @param A
     *            first matrix
     * @param D
     *            second matrix
     * @return second matrix
     *
     */
    public static <T extends Vector> T minusInplace(T A, Vector D) {
        for (int i = 0; i < A.size(); i++) {
            A.put(i, A.get(i) - D.get(i));
        }
        return A;
    }

    /**
     * Subtract a vector from another vector <code>A = A + D</code>
     *
     * @param A
     *            first matrix
     * @param D
     *            second matrix
     * @return second matrix
     *
     */
    public static <T extends Vector> T plusInplace(T A, Vector D) {
        for (int i = 0; i < A.size(); i++) {
            A.put(i, A.get(i) + D.get(i));
        }
        return A;
    }

    /**
     * Add two matrices, storing the results in the second:
     * <code>A = D + A</code>
     *
     * @param D
     *            first matrix
     * @param A
     *            second matrix
     * @return second matrix
     *
     */
    public static <T extends KuhnMatrix> T plusInplace(DiagonalMatrix D, T A) {
        double[] Dval = D.getVals();
        for (int i = 0; i < Dval.length; i++) {
            A.put(i, i, A.get(i, i) + Dval[i]);
        }
        return A;
    }

    /**
     * Compute A^T . B
     *
     * @param A
     * @param B
     * @return A^T . B
     */
    public static KuhnMatrix transposeDotProduct(KuhnMatrix A, KuhnMatrix B) {
        int mA = A.columnCount();
        int nB = B.columnCount();
        KuhnMatrix ret = A.newInstance(mA, nB);
        for (int i = 0; i < mA; i++) {
            Vector column = A.column(i);
            for (int j = 0; j < nB; j++) {
                double dot = column.dot(B.column(j));
                if (Math.abs(dot) > EPS) {
                    ret.put(i, j, dot);
                }
            }
        }
        return ret;
    }

    /**
     * Compute Y = A . B^T
     *
     * @param A
     * @param B
     * @return Y
     */
    public static KuhnMatrix dotProductTranspose(KuhnMatrix A, KuhnMatrix B) {
        KuhnMatrix ret = A.newInstance(A.rowCount(), B.rowCount());
        return dotProductTranspose(A, B, ret);
    }

    /**
     * Perform: A.T.dot(B.T) without performing the transpose. This is fine for
     * dense matrices but is very inefficient for sparse matrices, consider
     * performing the transpose manually.
     *
     * @param A
     * @param B
     * @return A.T.dot(B.T)
     */
    public static KuhnMatrix dotProductTransposeTranspose(KuhnMatrix A, KuhnMatrix B) {
        int mA = A.columnCount();
        int nB = B.rowCount();
        KuhnMatrix ret = A.newInstance(mA, nB);
        for (int i = 0; i < mA; i++) {
            Vector column = A.column(i);
            for (int j = 0; j < nB; j++) {
                double dot = column.dot(B.row(j));
                if (Math.abs(dot) > EPS) {
                    ret.put(i, j, dot);
                }
            }
        }
        return ret;
    }

    /**
     * Y = A . Bt
     *
     * @param A
     * @param B
     * @param Y
     * @return Y
     */
    public static <T extends KuhnMatrix> T dotProductTranspose(KuhnMatrix A, KuhnMatrix B,
            T Y)
    {
        if (A.columnCount() != B.columnCount()) {
            throw new RuntimeException(
                    String.format("Matrix size mismatch, A.cols == %d and B.T.cols == %d", A.columnCount(),
                            B.columnCount()));
        }
        int mA = A.rowCount();
        int nB = B.rowCount();
        for (int i = 0; i < mA; i++) {
            Vector arow = A.row(i);
            for (int j = 0; j < nB; j++) {
                Vector brow = B.row(j);
                double dot = arow.dot(brow);
                if (Math.abs(dot) > EPS) {
                    Y.put(i, j, dot);
                }
            }
        }
        return Y;
    }

    /**
     * A = A . s
     *
     * @param A
     * @param s
     * @return A
     */
    public static <T extends KuhnMatrix> T scaleInplace(T A, double s) {
        for (Vector row : A.rows()) {
            row.timesEquals(s);
        }
        return A;
    }

    /**
     * @param laplacian
     * @return returns a dense Jama matrix
     */
    public static JamaMatrix toJama(KuhnMatrix laplacian) {
        double[][] asArray = null;
        if (laplacian instanceof KuhnDenseMatrix) {
            asArray = ((KuhnDenseMatrix) laplacian).unwrap();
        } else {
            asArray = laplacian.asArray();
        }
        JamaMatrix ret = new JamaMatrix(asArray, laplacian.rowCount(),
                laplacian.columnCount());
        return ret;
    }

    /**
     * @param vector
     * @return the vector as a column in a matrix
     */
    public static JamaMatrix toColJama(Vector vector) {
        double[] vec = new double[vector.size()];
        vector.storeOn(vec, 0);
        JamaMatrix ret = new JamaMatrix(vec.length, 1);
        for (int i = 0; i < vec.length; i++) {
            ret.set(i, 0, vec[i]);
        }
        return ret;
    }

    /**
     * @param vector
     * @return the vector as a row in a matrix
     */
    public static JamaMatrix toRowJama(Vector vector) {
        double[] vec = new double[vector.size()];
        vector.storeOn(vec, 0);
        JamaMatrix ret = new JamaMatrix(1, vec.length);
        for (int i = 0; i < vec.length; i++) {
            ret.set(0, i, vec[i]);
        }
        return ret;
    }

    /**
     * @param sol
     * @return Dense matrix from a {@link gov.nist.math.jama.JamaMatrix}
     */
    public static KuhnMatrix fromJama(JamaMatrix sol) {
        return new KuhnDenseMatrix(sol.getArray());
    }

    /**
     * @param sol
     * @return Dense matrix from a {@link KuhnMatrix}
     */
    public static Matrix toMTJ(KuhnMatrix sol) {
        Matrix mat;
        if (sol instanceof SparseMatrix) {
            FlexCompRowMatrix fmat = new FlexCompRowMatrix(
                    sol.rowCount(), sol.columnCount());
            int i = 0;
            for (Vector vec : sol.rows()) {

                SparseVector x = new SparseVector(
                        vec.size(), vec.used());
                for (Entry ve : vec.entries()) {
                    x.set(ve.index, ve.value);
                }
                fmat.setRow(i, x);
                i++;
            }
            mat = fmat;
        } else {
            mat = new DenseMatrix(sol.rowCount(),
                    sol.columnCount());
            for (int i = 0; i < sol.rowCount(); i++) {
                for (int j = 0; j < sol.columnCount(); j++) {
                    mat.set(i, j, sol.get(i, j));
                }
            }
        }
        return mat;
    }

    /**
     * Extract the submatrix of the same type of mat
     *
     * @param mat
     * @param rowstart
     * @param rowend
     * @param colstart
     * @param colend
     * @return new instance
     */
    public static <T extends KuhnMatrix> T subMatrix(T mat, int rowstart,
            int rowend, int colstart, int colend)
    {
        @SuppressWarnings("unchecked")
        T ret = (T) mat.newInstance(rowend - rowstart, colend - colstart);

        for (int i = 0; i < ret.rowCount(); i++) {
            Vector row = mat.row(i + rowstart);
            for (Entry ent : row.entries()) {
                if (ent.index >= colstart && ent.index < colend) {
                    ret.put(i, ent.index - colstart, ent.value);
                }
            }
        }

        return ret;
    }

    /**
     * Calculate all 3, used by {@link #min(KuhnMatrix)}, {@link #max(KuhnMatrix)} and
     * {@link #mean(KuhnMatrix)}
     *
     * @param mat
     * @return the min, max and mean of the provided matrix
     */
    public static double[] minmaxmean(KuhnMatrix mat) {
        double min = Double.MAX_VALUE, max = -Double.MAX_VALUE, mean = 0;
        double size = mat.rowCount() * mat.columnCount();
        for (Vector v : mat.rows()) {
            for (Entry ent : v.entries()) {
                min = Math.min(min, ent.value);
                max = Math.max(max, ent.value);
                mean += ent.value / size;
            }
        }
        return new double[] { min, max, mean };
    }

    /**
     * Uses the first value returned by {@link #minmaxmean(KuhnMatrix)}
     *
     * @param mat
     * @return the min
     */
    public static double min(KuhnMatrix mat) {
        return minmaxmean(mat)[0];
    }

    /**
     * Uses the second value returned by {@link #minmaxmean(KuhnMatrix)}
     *
     * @param mat
     * @return the max
     */
    public static double max(KuhnMatrix mat) {
        return minmaxmean(mat)[1];
    }

    /**
     * Uses the third value returned by {@link #minmaxmean(KuhnMatrix)}
     *
     * @param mat
     * @return the mean
     */
    public static double mean(KuhnMatrix mat) {
        return minmaxmean(mat)[2];
    }

    /**
     * @param l
     * @param v
     * @return performs l - v returning a matrix of type T
     */
    public static <T extends KuhnMatrix> T minus(T l, double v) {
        @SuppressWarnings("unchecked")
        T ret = (T) l.newInstance(l.rowCount(), l.columnCount());
        int r = 0;
        for (Vector vec : l.rows()) {
            for (Entry ent : vec.entries()) {
                ret.put(r, ent.index, ent.value - v);
            }
            r++;
        }
        return ret;
    }

    /**
     * @param l
     * @param v
     * @return performs l - v returning a matrix of type T
     */
    public static Vector minus(Vector l, Vector v) {
        Vector ret = DenseVector.dense(l.size());
        for (int i = 0; i < l.size(); i++) {
            ret.put(i, l.get(i) - v.get(i));
        }
        return ret;
    }

    /**
     * @param v
     * @param l
     * @return performs v - l returning a matrix of type T
     */
    public static <T extends KuhnMatrix> T minus(double v, T l) {
        @SuppressWarnings("unchecked")
        T ret = (T) l.newInstance(l.rowCount(), l.columnCount());
        for (int i = 0; i < l.rowCount(); i++) {
            for (int j = 0; j < l.columnCount(); j++) {
                ret.put(i, j, v - l.get(i, j));
            }
        }
        return ret;
    }

    /**
     * Sum the diagonal of the given matrix
     *
     * @param d
     *            the matrix
     * @return the sum along the diagonal
     */
    public static double sum(DiagonalMatrix d) {
        double sum = 0.0;
        for (double v : d.getVals()) {
            sum += v;
        }
        return sum;
    }

    /**
     * Set a submatrix of a larger matrix
     *
     * @param to
     *            the matrix to write into
     * @param row
     *            the row to start inserting from
     * @param col
     *            the column to start inserting from
     * @param from
     *            the matrix to insert
     */
    public static void setSubMatrix(KuhnMatrix to, int row, int col, KuhnMatrix from) {
        for (int i = row; i < row + from.rowCount(); i++) {
            for (int j = col; j < col + from.columnCount(); j++) {
                to.put(i, j, from.get(i - row, j - col));
            }
        }
    }

    /**
     * Transpose a matrix, returning a new matrix.
     *
     * @param mat
     *            the matrix to transpose
     * @return the transposed matrix
     */
    public static <T extends KuhnMatrix> T transpose(T mat) {
        @SuppressWarnings("unchecked")
        T ret = (T) mat.newInstance(mat.columnCount(), mat.rowCount());
        for (int i = 0; i < mat.rowCount(); i++) {
            Vector v = mat.row(i);
            for (Entry ent : v.entries()) {
                ret.put(ent.index, i, ent.value);
            }
        }
        return ret;
    }

    /**
     * @param A
     * @param B
     * @return A = MAX(A,B)
     */
    public static SparseMatrix maxInplace(SparseMatrix A, SparseMatrix B) {
        for (int i = 0; i < A.rowCount(); i++) {
            Vector arow = A.row(i);
            Vector brow = B.row(i);
            for (Entry br : brow.entries()) {
                if (arow.get(br.index) < br.value) {
                    A.put(i, br.index, br.value);
                }
            }
        }
        return A;
    }

    /**
     * @param A
     * @param B
     * @return A = MIN(A,B)
     */
    public static SparseMatrix minInplace(SparseMatrix A, SparseMatrix B) {
        for (int i = 0; i < A.rowCount(); i++) {
            Vector arow = A.row(i);
            Vector brow = B.row(i);
            for (Entry br : brow.entries()) {
                if (arow.get(br.index) > br.value) {
                    A.put(i, br.index, br.value);
                }
            }
        }
        return A;
    }

    /**
     * @param A
     * @param B
     * @return A = A.*B
     */
    public static SparseMatrix timesInplace(SparseMatrix A, SparseMatrix B) {
        for (int i = 0; i < A.rowCount(); i++) {
            Vector arow = A.row(i);
            Vector brow = B.row(i);
            for (Entry br : brow.entries()) {
                A.put(i, br.index, br.value * arow.get(br.index));
            }
            for (Entry ar : arow.entries()) {
                 // All items in A not in B must be set to 0
                if (brow.get(ar.index) == 0.0) {
                    A.put(i, ar.index, 0.0);
                }
            }
        }
        return A;
    }

    /**
     * @param A
     * @param B
     * @return A = A.*B
     */
    public static KuhnMatrix timesInplace(KuhnMatrix A, KuhnMatrix B) {
        for (int i = 0; i < A.rowCount(); i++) {
            Vector arow = A.row(i);
            Vector brow = B.row(i);
            for (Entry br : brow.entries()) {
                A.put(i, br.index, br.value * arow.get(br.index));
            }
            for (Entry ar : arow.entries()) {
                // All items in A not in B must be set to 0
                if (brow.get(ar.index) == 0.0) {
                    A.put(i, ar.index, 0.0);
                }
            }
        }
        return A;
    }

    /**
     * Copy a matrix
     *
     * @param sparseMatrix
     *            the matrix to copy
     * @return the copy
     */
    public static <T extends KuhnMatrix> T copy(T sparseMatrix) {
        @SuppressWarnings("unchecked")
        T t = (T) sparseMatrix.newInstance();

        for (int r = 0; r < sparseMatrix.rowCount(); r++) {
            Vector row = sparseMatrix.row(r);
            for (Entry ent : row.entries()) {
                t.put(r, ent.index, ent.value);
            }
        }
        return t;
    }

    /**
     * Set values below the given threshold to zero in the output matrix.
     *
     * @param data
     *            the input matrix
     * @param thresh
     *            the threshold
     * @return a new matrix with values in the input matrix set to zero.
     */
    public static SparseMatrix threshold(SparseMatrix data, double thresh) {
        SparseMatrix newdata = (SparseMatrix) data.newInstance();
        for (int r = 0; r < data.rowCount(); r++) {
            Vector vec = data.row(r);
            for (Entry ent : vec.entries()) {
                if (ent.value > thresh) {
                    newdata.put(r, ent.index, 1);
                }
            }
        }
        return newdata;
    }

    /**
     *
     * @param to
     *            add items to this
     * @param startindex
     *            starting index in to
     * @param from
     *            add items from this
     */
    public static void setSubVector(Vector to, int startindex, Vector from) {
        if (to instanceof DenseVector && from instanceof DenseVector) {
            double[] tod = ((DenseVector) to).unwrap();
            double[] fromd = ((DenseVector) from).unwrap();
            System.arraycopy(fromd, 0, tod, startindex, fromd.length);
            return;
        }
        for (int i = 0; i < from.size(); i++) {
            to.put(i + startindex, from.get(i));
        }
    }

    /**
     * Starting from a given column of a row, set the values of a matrix to the
     * values of v
     *
     * @param to
     * @param row
     * @param col
     * @param v
     */
    public static void setSubMatrixRow(KuhnMatrix to, int row, int col, Vector v) {
        for (int i = col, j = 0; i < col + v.size(); i++, j++) {
            to.put(row, i, v.get(j));
        }
    }

    /**
     * Starting from a given row of a column, set the values of a matrix to the
     * values of v
     *
     * @param to
     * @param row
     * @param col
     * @param v
     */
    public static void setSubMatrixCol(KuhnMatrix to, int row, int col, Vector v) {
        for (int i = row, j = 0; i < row + v.size(); i++, j++) {
            to.put(i, col, v.get(j));
        }
    }

    /**
     * @param v
     *            the value vector
     * @param d
     *            the check value
     * @return for each item in the vector, returns 1 if the value is less than
     *         the check value
     */
    public static Vector lessThan(Vector v, double d) {
        Vector out = new KuhnSparseVector(v.size(), 1);
        for (int i = 0; i < v.size(); i++) {
            if (v.get(i) < d) {
                out.put(i, 1);
            }
        }
        return out;
    }

    /**
     * @param v
     * @return any values of the vector are non-zero
     */
    public static boolean any(Vector v) {
        for (int i = 0; i < v.size(); i++) {
            if (v.get(i) != 0.0) {
                return true;
            }
        }
        return false;
    }

    /**
     * @param m
     * @param col
     * @return m with the added column
     */
    public static <T extends KuhnMatrix> T appendColumn(T m, Vector col) {
        @SuppressWarnings("unchecked")
        T ret = (T) m.newInstance(m.rowCount(), m.columnCount() + 1);
        setSubMatrixCol(ret, 0, m.columnCount(), col);
        return ret;
    }

    /**
     * @param m
     * @param row
     * @return m with the added column
     */
    public static <T extends KuhnMatrix> T appendRow(T m, Vector row) {
        @SuppressWarnings("unchecked")
        T ret = (T) m.newInstance(m.rowCount() + 1, m.columnCount());
        setSubMatrixRow(ret, m.rowCount(), 0, row);
        return ret;
    }

    /**
     * Compute the dot product X.W
     *
     * @param X
     * @param W
     * @return dot product
     */
    public static KuhnMatrix dotProduct(KuhnMatrix X, KuhnMatrix W) {
        KuhnMatrix ret = X.newInstance(X.rowCount(), W.columnCount());
        for (int j = 0; j < ret.columnCount(); j++) {
            Vector column = W.column(j);
            for (int i = 0; i < ret.rowCount(); i++) {
                ret.put(i, j, X.row(i).dot(column));
            }
        }
        return ret;
    }

    /**
     * Compute the 2-norm (Euclidean norm) of the vector
     *
     * @param row
     *            the vector
     * @return the Euclidean norm
     */
    public static double norm2(Vector row) {
        double norm = 0.0;
        for (Entry e : row.entries()) {
            norm += (e.value * e.value);
        }
        return Math.sqrt(norm);
    }

    /**
     * Subtract matrices A-B
     *
     * @param A
     * @param B
     * @return A-B
     */
    public static KuhnMatrix minus(KuhnMatrix A, KuhnMatrix B) {
        KuhnMatrix ret = copy(A);
        minusInplace(ret, B);
        return ret;
    }

    /**
     * Compute the Frobenius norm
     *
     * @param A
     *            the matrix
     * @return the F norm
     */
    public static double normF(KuhnMatrix A) {
        double scale = 0.0, ssq = 1.0;
        for (Vector v : A.rows()) {
            for (Entry e : v.entries()) {
                double Aval = e.value;
                if (Aval != 0.0) {
                    double absxi = Math.abs(Aval);
                    if (scale < absxi) {
                        ssq = 1 + ssq * Math.pow(scale / absxi, 2.0);
                        scale = absxi;
                    } else {
                        ssq = ssq + Math.pow(absxi / scale, 2.0);
                    }
                }
            }
        }
        return scale * Math.sqrt(ssq);
    }

    /**
     * Stack matrices vertically
     *
     * @param matrices
     *            matrices to stack
     * @return matrix created from the stacking
     */
    public static KuhnMatrix vstack(KuhnMatrix... matrices) {
        int nrows = 0;
        int ncols = 0;
        for (KuhnMatrix matrix : matrices) {
            nrows += matrix.rowCount();
            ncols = matrix.columnCount();
        }
        KuhnMatrix ret = matrices[0].newInstance(nrows, ncols);
        int currentRow = 0;
        for (KuhnMatrix matrix : matrices) {
            setSubMatrix(ret, currentRow, 0, matrix);
            currentRow += matrix.rowCount();
        }
        return ret;
    }
}
