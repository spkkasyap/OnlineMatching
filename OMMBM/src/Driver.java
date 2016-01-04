/** This is the driver class from where our linkage and execution
 * to every other class or method in the project is connected. 
 * 
 * @author prathyush
 */

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class Driver {
	
	public static void generateSingleRun(int numNodes, String dataSource, double constant) {
		DateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
		Calendar cal = Calendar.getInstance();
		System.out.println("Start: " + dateFormat.format(cal.getTime()));
		System.out.println("Number of nodes in each set : "+numNodes);
		Matching m = new Matching();
		m.generateCostMatrix(numNodes);
		//m.printCostMatrix();
		System.out.println("The cost of matching produced by offline matching is : "+m.computeOfflineMatching(numNodes));
		
		//Calculating the cost of matching produced by Hungarian Algorithm
		double hungarian_cost = m.verifyHungarian();
		System.out.println("The cost of matching produced by Hungarian Algorithm is : "+hungarian_cost+"\n");
		
		cal = Calendar.getInstance();
		System.out.println("End: " + dateFormat.format(cal.getTime()));
}
	public static void main(String args[])
	{
		generateSingleRun(Integer.parseInt(args[0]), null, 0);
	}
}
