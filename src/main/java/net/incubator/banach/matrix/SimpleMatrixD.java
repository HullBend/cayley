package net.incubator.banach.matrix;

import java.util.Arrays;

public class SimpleMatrixD extends MatrixDBase implements MatrixD {

	public SimpleMatrixD(int rows, int cols) {
		super(rows, cols, new double[rows * cols], false);
	}

	public SimpleMatrixD(int rows, int cols, double initialValue) {
		super(rows, cols, new double[rows * cols], false);
		Arrays.fill(a, initialValue);
	}

	@Override
	public MatrixD multAdd(double alpha, MatrixD B, MatrixD C) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}

	@Override
	public MatrixD transABmultAdd(double alpha, MatrixD B, MatrixD C) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}

	@Override
	public MatrixD transAmultAdd(double alpha, MatrixD B, MatrixD C) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}

	@Override
	public MatrixD transBmultAdd(double alpha, MatrixD B, MatrixD C) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}
}
