package MDS;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class CityMDS {

	/**
	 * Reads the distances from file
	 * 
	 * @param path
	 * @return
	 * @throws IOException
	 */
	public static double[][] readDistances(String path) throws IOException {
		byte[] encoded = Files.readAllBytes(Paths.get(path));
		String datastr = new String(encoded, StandardCharsets.UTF_8);
		String[] rawData = datastr.split("\\s+");

		int size = (int) Math.sqrt(rawData.length) - 1;
		double[][] erg = new double[size][size];

		for (int i = 1; i <= size; i++) {
			for (int j = 1; j <= size; j++) {
				erg[i - 1][j - 1] = Double.parseDouble(rawData[i * (size + 1)
						+ j]);
			}
		}

		return erg;
	}

	/**
	 * Reads the name of the points from file
	 * 
	 * @param path
	 * @return
	 * @throws IOException
	 */
	public static String[] readCityNames(String path) throws IOException {
		byte[] encoded = Files.readAllBytes(Paths.get(path));
		String datastr = new String(encoded, StandardCharsets.UTF_8);
		String[] rawData = datastr.split("\\s+");

		int size = (int) Math.sqrt(rawData.length) - 1;
		String[] erg = new String[size];

		for (int i = 1; i <= size; i++) {
			erg[i - 1] = rawData[i];
		}

		return erg;
	}

	/**
	 * Visualization of the result
	 * 
	 * @param canvas
	 *            The canvas
	 * @param X
	 *            The coordinate matrix
	 * @param tags
	 *            The names of the points
	 */
	public static void showPoints(Draw canvas, RealMatrix X, String[] tags) {
		double xmin = Double.POSITIVE_INFINITY, xmax = Double.NEGATIVE_INFINITY, ymin = Double.POSITIVE_INFINITY, ymax = Double.NEGATIVE_INFINITY;

		// get canvas size
		for (int i = 0; i < X.getRowDimension(); i++) {
			double x = X.getEntry(i, 0);
			double y = X.getEntry(i, 1);

			if (x < xmin)
				xmin = x;
			if (x > xmax)
				xmax = x;
			if (y < ymin)
				ymin = y;
			if (y > ymax)
				ymax = y;
		}

		// System.out.println(xmin + " " + xmax);
		// System.out.println(ymin + " " + ymax);

		canvas.setXscale(xmin, xmax);
		canvas.setYscale(ymin, ymax);
		canvas.setPenRadius(0.01);

		for (int i = 0; i < X.getRowDimension(); i++) {
			canvas.point(X.getEntry(i, 0), X.getEntry(i, 1));
			canvas.text(X.getEntry(i, 0), X.getEntry(i, 1), tags[i]);
		}
	}

	/**
	 * Generates the squared distance matrix from the distance matrix
	 * 
	 * @param data
	 *            The distance matrix
	 * @return The squared distance matrix
	 */
	public static RealMatrix getSquaredDistanceMatrix(double[][] data) {
		// TODO
		double[][] value = new double[data.length][data[0].length];
		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data[i].length; j++) {
				value[i][j] = data[i][j] * data[i][j];
			}
		}
		RealMatrix valueMatrix = MatrixUtils.createRealMatrix(value);
		return valueMatrix;
	}

	/**
	 * Generates the scalar product matrix B
	 * 
	 * @param d2Matrix
	 *            The squared distance matrix
	 * @return The scalar product matrix
	 */
	public static RealMatrix getScalarProductMatrix(RealMatrix d2Matrix) {
		// TODO

		// Berechnung von J
		// J = I - 1/n 1*1^T
		RealMatrix EINHEITSMATRIX = MatrixUtils
				.createRealIdentityMatrix(d2Matrix.getRowDimension());

		double[][] einsen = new double[d2Matrix.getRowDimension()][d2Matrix
				.getColumnDimension()];
		for (int i = 0; i < einsen.length; i++) {
			for (int j = 0; j < einsen[i].length; j++) {
				einsen[i][j] = 1/(double)d2Matrix.getRowDimension();
			}
		}
		// 1/n 1*1^T
		RealMatrix EINS = MatrixUtils.createRealMatrix(einsen);

		RealMatrix J = EINHEITSMATRIX.subtract(EINS);

		// -1/2 j d^2 j
		RealMatrix B = (J.multiply(d2Matrix.multiply(J))).scalarMultiply(-0.5);
	//	 System.out.println(matrixToString(B));

		return B;
	}

	public static String matrixToString(RealMatrix m) {
		String s = "";
		int row = m.getRowDimension();
		int column = m.getColumnDimension();

		for (int i = 0; i < row; i++) {
			s += "[";
			for (int j = 0; j < column; j++) {
				s += m.getEntry(i, j) + " ";
			}
			s += "]\n";
		}
		return s;
	}

	/**
	 * Calculates the coordinate matrix X
	 * 
	 * @param B
	 *            The scalar product matrix
	 * @param eps
	 *            Thresshold for eigenvalues. Eigenvalues smaller then eps are
	 *            considered as zero or negative
	 * @return the coordinate matrix
	 */
	public static RealMatrix getCoordinates(RealMatrix B, double eps) {
		// TODO
		EigenDecomposition eigen = new EigenDecomposition(B);
		// ermittele !0 Eigenwerte
		int counter =0;
		int[] notNullEigenValues = new int[B.getColumnDimension()];
		for (int i = 0; i <  B.getColumnDimension(); i++) {
			if ( Math.abs(eigen.getRealEigenvalue(i)) > eps )
				{
				counter++;
				notNullEigenValues[counter-1]= i;
				}
			
		//	System.out.println( ( eigen.getRealEigenvalue(i)));
		}
//		for (int i = 0; i < notNullEigenValues.length-1; i++) {
//			System.out.println(notNullEigenValues[i]);
//		}
//		System.out.println(counter);

		
		

		RealMatrix LAMBDA = MatrixUtils.createRealIdentityMatrix(counter);
		for (int i = 0; i < counter; i++) {
		
			LAMBDA.setEntry(i, i, Math.pow(eigen.getRealEigenvalue(notNullEigenValues[i]), 0.5));

		}
		// System.out.println(matrixToString(LAMBDA));
		 

		 RealMatrix Q = MatrixUtils.createRealMatrix(B.getRowDimension(), counter);
		 for (int i = 0; i < counter; i++) {
			 
			Q.setColumnVector(i, eigen.getEigenvector(notNullEigenValues[i]));
		
			 
		}
		 //System.out.println(matrixToString(Q));

		RealMatrix X = Q.multiply(LAMBDA);
		 System.out.println(matrixToString(X));

		return X;
	}

	public static void main(String[] args) throws IOException {
		// Read the measured data from file and the name of the points
	//	double[][] data = readDistances("res/staedtedistanzenluft.txt");
	//	String[] tags = readCityNames("res/staedtedistanzenluft.txt");

		 double[][] data = readDistances("res/mdstestdata.txt");
		 String[] tags = readCityNames("res/mdstestdata.txt");

		RealMatrix dsquare = getSquaredDistanceMatrix(data);
		RealMatrix B = getScalarProductMatrix(dsquare);
		RealMatrix X = getCoordinates(B, 1.0E-6);
	//	System.out.println(X);

		// Visualize data
		Draw canvas = new Draw();
		showPoints(canvas, X, tags);
	}
}
