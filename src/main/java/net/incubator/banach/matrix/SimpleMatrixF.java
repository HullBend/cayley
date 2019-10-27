package net.incubator.banach.matrix;

import java.util.Arrays;

public class SimpleMatrixF extends MatrixFBase implements MatrixF {

	public SimpleMatrixF(int rows, int cols) {
		super(rows, cols, new float[rows * cols], false);
	}

	public SimpleMatrixF(int rows, int cols, float initialValue) {
		super(rows, cols, new float[rows * cols], false);
		Arrays.fill(a, initialValue);
	}

	@Override
	public MatrixF multAdd(float alpha, MatrixF B, MatrixF C) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}

	@Override
	public MatrixF transABmultAdd(float alpha, MatrixF B, MatrixF C) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}

	@Override
	public MatrixF transAmultAdd(float alpha, MatrixF B, MatrixF C) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}

	@Override
	public MatrixF transBmultAdd(float alpha, MatrixF B, MatrixF C) {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}
}
