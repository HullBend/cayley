package ch.akuhn.matrix;

import java.util.Arrays;

/**
 * Matrix where <code>a<sub>ij</sub> = a<sub>ji</sub></code> for all elements.
 * <p>
 * 
 * @author Adrian Kuhn
 */
public class SymmetricMatrix extends KuhnDenseMatrix {

    /**
     * Construct with given size
     * 
     * @param size
     */
    public SymmetricMatrix(int size) {
        super(size, size);
    }

    /**
     * Construct with given values, which must be jagged and represent the lower
     * triangular values
     * 
     * @param values
     */
    public SymmetricMatrix(double[][] values) {
        super(values);
    }

    @Override
    protected void assertInvariant() throws IllegalArgumentException {
        for (int n = 0; n < values.length; n++) {
            if (values[n].length != (n + 1)) {
                throw new IllegalArgumentException();
            }
        }
    }

    @Override
    protected double[][] makeValues(int rows, int columns) {
        if (rows != columns) {
            throw new IllegalArgumentException("rows != columns");
        }
        double[][] values = new double[rows][];
        for (int n = 0; n < values.length; n++) {
            values[n] = new double[n + 1];
        }
        return values;
    }

    @Override
    public int columnCount() {
        return rowCount();
    }

    @Override
    public double get(int row, int column) {
        return row > column ? values[row][column] : values[column][row];
    }

    @Override
    public double put(int row, int column, double value) {
        return row > column ? (values[row][column] = value) : (values[column][row] = value);
    }

    @Override
    public int rowCount() {
        return values.length;
    }

    /**
     * Create from a square matrix
     * 
     * @param square
     * @return the matrix
     */
    public static KuhnDenseMatrix fromSquare(double[][] square) {
        double[][] jagged = new double[square.length][];
        for (int i = 0; i < jagged.length; i++) {
            if (square[i].length != square.length) {
                throw new IllegalArgumentException("square[i].length != square.length : i = " + i);
            }
            jagged[i] = Arrays.copyOf(square[i], i + 1);
        }
        return new SymmetricMatrix(jagged);
    }

    /**
     * Create from jagged low triangular values
     * 
     * @param values
     * @return the matrix
     */
    public static KuhnDenseMatrix fromJagged(double[][] values) {
        return new SymmetricMatrix(values);
    }

    @Override
    public double[][] unwrap() {
        return values;
    }

    @Override
    public double[] rowwiseMean() {
        double[] mean = new double[rowCount()];
        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < i; j++) {
                mean[i] += values[i][j];
                mean[j] += values[i][j];
            }
        }
        for (int n = 0; n < mean.length; n++) {
            mean[n] /= mean.length;
        }
        return mean;
    }

    @Override
    public Vector mult(Vector v) {
        if (v.size() != values.length) {
            throw new IllegalArgumentException("Vector.size() != values.length : " + v.size());
        }
        double[] mult = new double[v.size()];
        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < i; j++) {
                mult[i] += values[i][j] * v.get(j);
                mult[j] += values[i][j] * v.get(i);
            }
        }
        return Vector.wrap(mult);
    }
}
