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

import java.util.ArrayList;
import java.util.List;

/**
 * {@code IndependentPair} represents a generic pair of objects of different
 * (independent) types.
 * 
 * @author Jonathon Hare
 * 
 * @param <A>
 *            the class of the first object in the pair
 * @param <B>
 *            the class of the second object in the pair
 */
public class IndependentPair<A, B> {

    private final A o1;
    private final B o2;

    /**
     * Constructs a Pair object with two objects obj1 and obj2
     * 
     * @param obj1
     *            first object in pair
     * @param obj2
     *            second object in pair
     */
    public IndependentPair(A obj1, B obj2)
    {
        this.o1 = obj1;
        this.o2 = obj2;
    }

    /**
     * @return first object in pair
     */
    public A firstObject()
    {
        return o1;
    }

    /**
     * @return second object in pair
     */
    public B secondObject()
    {
        return o2;
    }

    /**
     * @return first object in pair
     */
    public A getFirstObject()
    {
        return o1;
    }

    /**
     * @return second object in pair
     */
    public B getSecondObject()
    {
        return o2;
    }

    @Override
    public String toString() {
        return "[" + this.o1 + "," + this.o2 + "]";
    }

    @Override
    public boolean equals(Object thatObject) {
        if (!(thatObject instanceof IndependentPair)) {
            return false;
        }
        @SuppressWarnings("rawtypes")
        IndependentPair that = (IndependentPair) thatObject;
        return this.o1 == that.o1 && this.o2 == that.o2;
    }

    /**
     * Create a pair from the given objects.
     * 
     * @param <T>
     *            Type of first object.
     * @param <Q>
     *            Type of second object.
     * @param t
     *            The first object.
     * @param q
     *            The second object.
     * @return The pair.
     */
    public static <T, Q> IndependentPair<T, Q> pair(T t, Q q) {
        return new IndependentPair<T, Q>(t, q);
    }

    /**
     * Extract the first objects from a list of pairs.
     * 
     * @param <T>
     *            type of first object
     * @param <Q>
     *            type of second object
     * @param data
     *            the data
     * @return extracted first objects
     */
    public static <T, Q> List<T> getFirst(Iterable<? extends IndependentPair<T, Q>> data) {
        ArrayList<T> extracted = new ArrayList<T>();

        for (IndependentPair<T, Q> item : data) {
            extracted.add(item.o1);
        }

        return extracted;
    }

    /**
     * Extract the second objects from a list of pairs.
     * 
     * @param <T>
     *            type of first object
     * @param <Q>
     *            type of second object
     * @param data
     *            the data
     * @return extracted second objects
     */
    public static <T, Q> List<Q> getSecond(Iterable<? extends IndependentPair<T, Q>> data) {
        ArrayList<Q> extracted = new ArrayList<Q>();

        for (IndependentPair<T, Q> item : data) {
            extracted.add(item.o2);
        }

        return extracted;
    }

    /**
     * Get the function that returns the first object from the pair
     * 
     * @return the function that returns the first object from the pair
     */
    public static <T, Q> MathFunction<IndependentPair<T, Q>, T> getFirstFunction() {
        return new MathFunction<IndependentPair<T, Q>, T>() {
            @Override
            public T apply(IndependentPair<T, Q> in) {
                return in.o1;
            }
        };
    }

    /**
     * Get the function that returns the second object from the pair
     * 
     * @return the function that returns the second object from the pair
     */
    public static <T, Q> MathFunction<IndependentPair<T, Q>, Q> getSecondFunction() {
        return new MathFunction<IndependentPair<T, Q>, Q>() {
            @Override
            public Q apply(IndependentPair<T, Q> in) {
                return in.o2;
            }
        };
    }

    /**
     * Create a pair list from the given objects.
     * 
     * @param <T>
     *            Type of objects.
     * @param t
     *            The list of first objects.
     * @param q
     *            The list of second objects.
     * @return The list of pairs.
     */
    public static <T, Q> List<IndependentPair<T, Q>> pairList(List<T> t, List<Q> q) {
        ArrayList<IndependentPair<T, Q>> list = new ArrayList<IndependentPair<T, Q>>(t.size());

        for (int i = 0; i < t.size(); i++) {
            list.add(new IndependentPair<T, Q>(t.get(i), q.get(i)));
        }

        return list;
    }

    /**
     * Create a new {@link IndependentPair} from this one with the elements
     * swapped
     * 
     * @return the swapped pair
     */
    public IndependentPair<B, A> swap() {
        return new IndependentPair<B, A>(o2, o1);
    }

    /**
     * Swap the order of the pairs
     * 
     * @param data
     *            the input
     * @return the swapped data
     */
    public static <T, Q> List<IndependentPair<? extends Q, ? extends T>> swapList(
            List<? extends IndependentPair<? extends T, ? extends Q>> data)
    {
        ArrayList<IndependentPair<? extends Q, ? extends T>> list = new ArrayList<IndependentPair<? extends Q, ? extends T>>(
                data.size());

        for (int i = 0; i < data.size(); i++) {
            list.add(data.get(i).swap());
        }

        return list;
    }
}
