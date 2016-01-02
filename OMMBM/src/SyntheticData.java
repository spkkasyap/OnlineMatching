import java.util.ArrayList;
import java.util.Random;

/**
 * Generates random synthetic 1D and 2D data and the associated cost matrix.
 */
public class SyntheticData {
	
	public static int[][] coefficientMapping;

	/**
	 * Generates 1D synthetic data and creates a cost matrix. Picks a random
	 * integer between 0 and 10000 and adds it to setA (no repeats). Same thing
	 * for setB. Generates the absolute value distance between each node of setA
	 * and setB and creates a cost matrix.
	 * 
	 * @param numSetA
	 *            The number of taxis (nodes in set A)
	 * @return Cost matrix of the data
	 */
	public double[][] generateSynthetic1D(int numSetA) {
		
		coefficientMapping = new int[numSetA][numSetA];

		Random randSetA = new Random();
		Random randSetB = new Random();

		ArrayList<Integer> setA = new ArrayList<Integer>();
		ArrayList<Integer> setB = new ArrayList<Integer>();

		for (int i = 0; i < numSetA; i++) {
			while (true) {
				int test = randSetA.nextInt(10000);

				if (!setA.contains(test)) {
					setA.add(test);
				} else {
					continue;
				}
				break;
			}
			while (true) {
				int test2 = randSetB.nextInt(10000);

				if (!setB.contains(test2)) {
					setB.add(test2);
				} else {
					continue;
				}
				break;
			}
		}
		
		for (int i = 0; i < numSetA; i++) {
			for (int j = 0; j < numSetA; j++) {
				coefficientMapping[i][j] = setA.get(i) < setB.get(j) ? 3 : 6;
			}
		}

		double[][] costMatrix = new double[numSetA][numSetA];

		for (int i = 0; i < costMatrix.length; i++) {
			for (int j = 0; j < costMatrix[i].length; j++) {
				costMatrix[i][j] = Math.abs(setA.get(i) - setB.get(j));
			}
		}
		return costMatrix;
	}

}
