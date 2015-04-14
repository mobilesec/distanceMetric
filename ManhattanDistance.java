

/**
 * 
 * This class is an implimentation of Manhatten distance 
 * @version 1.0
 * @author Muhammad Muaaz
 * 
 */
public class ManhattanDistance {


	public static double calculateManHattenDistance(double[]array1, double[] array2) throws IllegalArgumentException{

		double sum = 0.0;
		if(array1.length == array2.length)
		{
			for (int i = 0; i <array1.length; i++)
			{
				sum += Math.abs(array1[i] - array2[i]);
			}
		}
		else
		{
			throw new IllegalArgumentException("Both arrays should have equal length");
		}
		return sum;

	}

}
