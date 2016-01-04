/** This is the driver class from where our linkage and execution
 * to every other class or method in the project is connected. 
 * 
 * @author prathyush
 */

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class Driver {
	private double cost_offline;
	private double cost_online;
	private double cost_hungarian;
	private double cost_greedy;
	private double cr_offbyon;
	
	/** A constructor for the driver class
	 * 
	 */
	public Driver(){
		this.cost_offline = 0.0;
		this.cost_online = 0.0;
		this.cost_hungarian = 0.0;
		this.cost_greedy = 0.0;
		this.cr_offbyon = 0.0;
	}
	
	/** A function to generate a single run of the algorithm 
	 * 
	 * @param numNodes
	 * @param dataSource
	 * @param constant
	 */
	public void generateSingleRun(int numNodes, String dataSource, double constant) {
		DateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
		Calendar cal = Calendar.getInstance();
		System.out.println("Start: " + dateFormat.format(cal.getTime()));
		System.out.println("Number of nodes in each set : "+numNodes);
		Matching m = new Matching();
		m.generateCostMatrix(numNodes);
		
		//m.printCostMatrix();
		this.cost_offline = m.computeOfflineMatching(numNodes);
		this.cost_online = m.computeOnlineMatching(numNodes, m.permuteDestinations(numNodes));
		this.cr_offbyon = (this.cost_online)/(this.cost_offline);
		
		System.out.println("The cost of matching produced by offline matching is : "+this.cost_offline);
		System.out.println("The cost of matching produced by online matching is :"+this.cost_online);
		System.out.println("Competitive Ratio: online/offline : "+this.cr_offbyon);
		
		//Calculating the cost of matching produced by Hungarian Algorithm
		this.cost_hungarian = m.verifyHungarian();
		System.out.println("The cost of matching produced by Hungarian Algorithm is : "+this.cost_hungarian+"\n");
		
		cal = Calendar.getInstance();
		System.out.println("End: " + dateFormat.format(cal.getTime()));
}
	public static void main(String args[])
	{
		Driver driver = new Driver();
		driver.generateSingleRun(Integer.parseInt(args[0]), null, 0);
	}
}
