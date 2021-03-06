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

import ch.akuhn.matrix.KuhnMatrix;
import gov.nist.math.jama.JamaMatrix;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * Methods for calculating the Moore-Penrose Pseudo-Inverse
 *
 * @author Jonathon Hare (jsh2@ecs.soton.ac.uk)
 */
public class PseudoInverse {

    /**
     * Compute the Moore-Penrose Pseudo-Inverse.
     *
     * @param matrix
     *            The matrix to invert.
     * @return the pseudo-inverse.
     */
    public static JamaMatrix pseudoInverse(JamaMatrix matrix) {
        final DenseMatrix mjtA = new DenseMatrix(matrix.getArray());
        SVD svd;

        try {
            svd = SVD.factorize(mjtA);
        } catch (NotConvergedException e) {
            throw new RuntimeException(e);
        }

        final JamaMatrix Sinv = new JamaMatrix(matrix.getColumnDimension(), matrix.getRowDimension());

        final double[] Sarr = svd.getS();
        for (int i = 0; i < svd.getS().length; i++) {
            if (Sarr[i] != 0.0) {
                Sinv.set(i, i, 1.0 / Sarr[i]);
            }
        }

        final JamaMatrix Vt = new JamaMatrix(svd.getVt().numRows(), svd.getVt().numColumns());
        for (int r = 0; r < svd.getVt().numRows(); r++) {
            for (int c = 0; c < svd.getVt().numColumns(); c++) {
                Vt.set(r, c, svd.getVt().get(r, c));
            }
        }

        final JamaMatrix U = new JamaMatrix(svd.getU().numRows(), svd.getU().numColumns());
        for (int r = 0; r < svd.getU().numRows(); r++) {
            for (int c = 0; c < svd.getU().numColumns(); c++) {
                U.set(r, c, svd.getU().get(r, c));
            }
        }

        final JamaMatrix pinv = Vt.transpose().times(Sinv).times(U.transpose());

        return pinv;
    }

    /**
     * Compute the lower-rank approximation of the Moore-Penrose Pseudo-Inverse.
     *
     * @param matrix
     *            The matrix to invert.
     * @param rank
     *            the desired rank.
     * @return the pseudo-inverse.
     */
    public static JamaMatrix pseudoInverse(JamaMatrix matrix, int rank) {
        return pseudoInverse(new JamaDenseMatrix(matrix), rank);
    }

    /**
     * Compute the lower-rank approximation of the Moore-Penrose Pseudo-Inverse.
     *
     * @param matrix
     *            The matrix to invert.
     * @param rank
     *            the desired rank.
     * @return the pseudo-inverse.
     */
    public static JamaMatrix pseudoInverse(KuhnMatrix matrix, int rank) {
        final ThinSingularValueDecomposition tsvd = new ThinSingularValueDecomposition(matrix, rank);

        final JamaMatrix Sinv = new JamaMatrix(tsvd.S.length, tsvd.S.length);
        for (int i = 0; i < tsvd.S.length; i++) {
            if (tsvd.S[i] != 0.0) {
                Sinv.set(i, i, 1.0 / tsvd.S[i]);
            }
        }

        final JamaMatrix pinv = tsvd.Vt.transpose().times(Sinv).times(tsvd.U.transpose());

        return pinv;
    }
}
