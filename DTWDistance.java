/*
 * DTW.java
 */

/**
 * This class implements the Dynamic Time Warping algorithm given two sequences. 
 * I have checked the out put of this class and it works exactly same as of the 
 * MATLAB dtw implementation. I have added one extra function to this class to 
 * get ACCUMULATED-DISTANCE. This accumulated distance is basically the distance
 * between two time series. The original implementation only returns the warping 
 * distance.
 * <pre>
 *   X = x1, x2,..., xi,..., xn
 *   Y = y1, y2,..., yj,..., ym
 *  </pre>
 *  
 * @author Cheol-Woo Jung (cjung@gatech.edu) Modified: Muhammad Muaaz 
 * @version	1.0	
 */
public class DTWDistance {

	protected double[] seq1;
	protected double[] seq2;
	protected int[][] warpingPath;

	protected int n;
	protected int m;
	protected int K;

	protected double warpingDistance;
	protected double accumulatedDist;

	/**
	 * Constructor
	 * 
	 * @param query		
	 * @param templete	
	 */
	public DTWDistance(double[] sample, double[] templete) {
		seq1 = sample;
		seq2 = templete;

		n = seq1.length;	
		m = seq2.length;
		K = 1;

		warpingPath = new int[n + m][2];	// max(n, m) <= K < n + m
		warpingDistance = 0.0;
		accumulatedDist = 0.0;

		this.compute();
	}

	public void compute() {
		double accumulatedDistance = 0.0;

		double[][] d = new double[n][m];	// local distances
		double[][] D = new double[n][m];	// global distances

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				d[i][j] = distanceBetween(seq1[i], seq2[j]);
			}
		}

		D[0][0] = d[0][0];

		for (int i = 1; i < n; i++) {
			D[i][0] = d[i][0] + D[i - 1][0];
		}

		for (int j = 1; j < m; j++) {
			D[0][j] = d[0][j] + D[0][j - 1];
		}

		for (int i = 1; i < n; i++) {
			for (int j = 1; j < m; j++) {
				accumulatedDistance = Math.min(Math.min(D[i-1][j], D[i-1][j-1]), D[i][j-1]);
				accumulatedDistance += d[i][j];
				D[i][j] = accumulatedDistance;
			}
		}
		accumulatedDistance = D[n - 1][m - 1];

		int i = n - 1;
		int j = m - 1;
		int minIndex = 1;

		warpingPath[K - 1][0] = i;
		warpingPath[K - 1][1] = j;

		while ((i + j) != 0) {
			if (i == 0) {
				j -= 1;
			} else if (j == 0) {
				i -= 1;
			} else {	// i != 0 && j != 0
				double[] array = { D[i - 1][j], D[i][j - 1], D[i - 1][j - 1] };
				minIndex = this.getIndexOfMinimum(array);

				if (minIndex == 0) {
					i -= 1;
				} else if (minIndex == 1) {
					j -= 1;
				} else if (minIndex == 2) {
					i -= 1;
					j -= 1;
				}
			} // end else
			K++;
			warpingPath[K - 1][0] = i;
			warpingPath[K - 1][1] = j;
		} // end while
		warpingDistance = accumulatedDistance / K;
		this.accumulatedDist = accumulatedDistance;
		this.reversePath(warpingPath);
	}

	/**
	 * Changes the order of the warping path (increasing order)
	 * 
	 * @param path	the warping path in reverse order
	 */
	protected void reversePath(int[][] path) {
		int[][] newPath = new int[K][2];
		for (int i = 0; i < K; i++) {
			for (int j = 0; j < 2; j++) {
				newPath[i][j] = path[K - i - 1][j];
			}
		}
		warpingPath = newPath;
	}
	/**
	 * Returns the warping distance
	 * 
	 * @return
	 */
	public double getWarpingDistance() {
		return warpingDistance;
	}

	/**
	 * Returns the warping distance
	 * 
	 * @return
	 */
	public double getAccumulatedDistance() {
		return accumulatedDist;
	}
	/**
	 * Computes a distance between two points
	 * 
	 * @param p1	the point 1
	 * @param p2	the point 2
	 * @return		the distance between two points
	 */
	protected double distanceBetween(double p1, double p2) {
		return (p1 - p2) * (p1 - p2);
	}

	/**
	 * Finds the index of the minimum element from the given array
	 * 
	 * @param array		the array containing numeric values
	 * @return				the min value among elements
	 */
	protected int getIndexOfMinimum(double[] array) {
		int index = 0;
		double val = array[0];

		for (int i = 1; i < array.length; i++) {
			if (array[i] < val) {
				val = array[i];
				index = i;
			}
		}
		return index;
	}

	/**
	 *	Returns a string that displays the warping distance and path
	 */
	public String toString() {
		String retVal = "Warping Distance: " + warpingDistance + "\n";
		retVal += "Warping Path: {";
		for (int i = 0; i < K; i++) {
			retVal += "(" + warpingPath[i][0] + ", " +warpingPath[i][1] + ")";
			retVal += (i == K - 1) ? "}" : ", ";

		}
		return retVal;
	}


/*	public double cyclicRotationMatric()throws IllegalArgumentException {
		double crmDist = 0.0;

		//check if length of reference and probe cycle is equal...

		if (seq1.length!=seq2.length) {
			throw new IllegalArgumentException("Length of Template and proble cycle must be equal for CRM");
		}
		double []distances = new double[seq1.length];
		for (int i = 0; i <= seq1.length; i++) {

			double [] rotatedSeq = new double [seq1.length];
			double firsVal = seq1[0];
			for (int j = 1; j <seq1.length; j++) {
				rotatedSeq[j-1]= seq1[j];
			}
			rotatedSeq[seq1.length] = firsVal;

			distances[i] = ManhattanDistance.calculateManHattenDistance(rotatedSeq, seq2);
		}

		crmDist = ArrayManupulation.minimum(distances);

		return crmDist;
	}
*/

	/**
	 * Tests functionality of this class
	 * 
	 * @param args	ignored
	 */
	public static void main(String[] args) {
		double[] n2 = {1.5, 3.9, 4.1, 3.3, 2.1, 2.45, 3.673};
		double[] n1 = {2.1, 2.45, 3.673, 4.32, 2.05, 1.93, 5.67, 6.01};
		DTWDistance dtw = new DTWDistance(n1, n2);
		System.out.println("Accumulated distance between the series is :  "+dtw.getAccumulatedDistance());
		System.out.println("Warping Distance is between the series is :  "+dtw.getWarpingDistance());
	}
}
