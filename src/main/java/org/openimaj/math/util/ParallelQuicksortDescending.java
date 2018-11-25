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
package org.openimaj.math.util;

/**
 * Sort parallel arrays. Arrays are sorted in-place. The first array determines
 * the order, and is sorted into descending order.
 * <p>
 * 
 * @author Jonathan Hare
 * @author Sina Samangooei
 */
public final class ParallelQuicksortDescending {

    /**
     * Sort parallel arrays. Arrays are sorted in-place. The first array
     * determines the order, and is sorted into descending order.
     * <p>
     * Implementation inspired by this stackoverflow page: <a href=
     * "http://stackoverflow.com/questions/951848/java-array-sort-quick-way-to-get-a-sorted-list-of-indices-of-an-array"
     * > http://stackoverflow.com/questions/951848/java-array-sort-quick-way-to-
     * get-a-sorted-list-of-indices-of-an-array </a>
     * 
     * @param main
     *            the values to use for determining the order
     * @param indices
     *            the second array
     */
    public static void sort(double[] main, int[] indices) {
        sort(main, indices, 0, indices.length - 1);
    }

    /**
     * Sort parallel arrays. Arrays are sorted in-place. The first array
     * determines the order, and is sorted into descending order.
     * <p>
     * Implementation inspired by this stackoverflow page: <a href=
     * "http://stackoverflow.com/questions/951848/java-array-sort-quick-way-to-get-a-sorted-list-of-indices-of-an-array"
     * > http://stackoverflow.com/questions/951848/java-array-sort-quick-way-to-
     * get-a-sorted-list-of-indices-of-an-array </a>
     * 
     * @param main
     *            the values to use for determining the order
     * @param indices
     *            the second array
     * @param left
     *            the starting index
     * @param right
     *            the ending index
     */
    public static void sort(double[] main, int[] indices, int left,
            int right)
    {
        if (right <= left) {
            return;
        }

        int i = partitionDesc(main, indices, left, right);

        sort(main, indices, left, i - 1);
        sort(main, indices, i + 1, right);
    }

    // partition a[left] to a[right], assumes left < right
    private static int partitionDesc(double[] a, int[] index, int left, int right) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (a[++i] > a[right])
                // find item on left to swap
                ; // a[right] acts as sentinel
            while (a[right] > a[--j])
                // find item on right to swap
                if (j == left) {
                    break; // don't go out-of-bounds
                }
            if (i >= j) {
                break; // check if pointers cross
            }
            exch(a, index, i, j); // swap two elements into place
        }
        exch(a, index, i, right); // swap with partition element
        return i;
    }

    // exchange a[i] and a[j]
    private static void exch(double[] a, int[] index, int i, int j) {
        double swap = a[i];
        a[i] = a[j];
        a[j] = swap;

        int b = index[i];
        index[i] = index[j];
        index[j] = b;
    }

    private ParallelQuicksortDescending() {
    }
}
