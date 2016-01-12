/** This is the driver class from where our linkage and execution
 * to every other class or method in the project is connected. 
 * 
 * @author prathyush
 */

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;

public class Driver {
	private double cost_offline;
	private double cost_online;
	private double cost_hungarian;
	private double cost_greedy;
	private double cr_onbyoff;
	private double cr_grbyoff;

	/** A constructor for the driver class
	 * 
	 */
	public Driver(){
		this.cost_offline = 0.0;
		this.cost_online = 0.0;
		this.cost_hungarian = 0.0;
		this.cost_greedy = 0.0;
		this.cr_onbyoff = 0.0;
		this.cr_grbyoff = 0.0;
	}

	/** A function to generate a single run of the algorithm 
	 * 
	 * @param numNodes
	 * @param dataSource
	 * @param constant
	 * @throws IOException 
	 */
	public void generateSingleRun(int numNodes, String dataSource, int constant) throws IOException {
		ArrayList<Integer> destinationIndices;
		BufferedWriter bw = new BufferedWriter(new FileWriter("Output.txt"));
		DateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
		Calendar cal = Calendar.getInstance();
		
		bw.write("Start: " + dateFormat.format(cal.getTime())+"\n");
		bw.write("Number of nodes in each set : "+numNodes+"\n");

		Matching m = new Matching();
		m.generateCostMatrix(numNodes);

		//m.printCostMatrix();
		destinationIndices = m.permuteDestinations(numNodes);
		
		//this.cost_offline = m.computeOfflineMatching(numNodes);
		this.cost_online = m.computeOnlineMatchingDW(numNodes, destinationIndices, constant);
		//this.cost_online = m.computeOnlineMatching(numNodes, destinationIndices);
		//this.cr_onbyoff = (this.cost_online)/(this.cost_offline);

		//bw.write("The cost of matching produced by offline matching is : "+this.cost_offline+"\n");
		bw.write("The cost of matching produced by online matching is :"+this.cost_online+"\n");
		//bw.write("Competitive Ratio: online/offline : "+this.cr_onbyoff+"\n");

		//Calculating the cost of matching produced by Hungarian Algorithm
		this.cost_hungarian = m.verifyHungarian();
		bw.write("The cost of matching produced by Hungarian Algorithm is : "+this.cost_hungarian+"\n");
		
		//Calculating the cost of online greedy matching
		this.cost_greedy = m.computeGreedyMatching(numNodes, destinationIndices);
		bw.write("The cost of matching produced by Online Greedy Algorithm is : "+this.cost_greedy+"\n");
		this.cr_grbyoff = (this.cost_greedy)/(this.cost_offline);
		bw.write("Competitive Ratio: greedy/offline : "+this.cr_grbyoff+"\n");
		
		cal = Calendar.getInstance();
		bw.write("End: " + dateFormat.format(cal.getTime())+"\n");
		
		bw.close();
		System.out.println("Check the output file");
	}


	public static void main(String args[]) throws NumberFormatException, IOException
	{
		Driver driver = new Driver();
		driver.generateSingleRun(Integer.parseInt(args[0]), null, Integer.parseInt(args[1]));
	}
}
