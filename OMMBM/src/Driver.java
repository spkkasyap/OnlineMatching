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
	private ArrayList<Integer> constants;

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
		this.constants = new ArrayList<Integer>();
		this.constants.add(1);
		this.constants.add(3);
		this.constants.add(5);
		this.constants.add(10);
		this.constants.add(100);
		this.constants.add(1000);
		this.constants.add(10000);
		this.constants.add(100000);
		this.constants.add(1000000);
	}

	/** A function to generate a single run of the algorithm 
	 * 
	 * @param numNodes
	 * @param dataSource
	 * @param constant
	 * @throws IOException 
	 */
	public String generateSingleRun(int numNodes, String dataSource) throws IOException {
		StringBuilder s = new StringBuilder();
		ArrayList<Integer> destinationIndices;
		DateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
		Calendar cal = Calendar.getInstance();

		System.out.println("Start: " + dateFormat.format(cal.getTime())+"\n");
		System.out.println("Number of nodes in each set : "+numNodes+"\n");


		//bw.write("Start: " + dateFormat.format(cal.getTime())+"\n");
		//bw.write("Number of nodes in each set : "+numNodes+"\n");

		Matching m = new Matching();
		m.generateCostMatrix(numNodes);

		//m.printCostMatrix();
		destinationIndices = m.permuteDestinations(numNodes);
		int run = 1;
		s.append(numNodes+", ");
		for(int k = 0; k < this.constants.size(); k++)
		{
			int t = this.constants.get(k);
			System.out.println("The constant t is "+t);
			s.append(m.computeOnlineMatchingDW(numNodes, destinationIndices, t));
			if(k < this.constants.size()-1)
				s.append(", ");
			//s.append("-----------------------different value of t---------------------\n");
			System.out.println("Run : "+run+" completed!");
			run++;
		}
		s.append("\n");

		//s.append("=================HUNGARIAN AND ONLINE GREEDY RESULTS=====================\n");
		this.cost_hungarian = m.verifyHungarian();
		//s.append("The cost of matching produced by Hungarian Algorithm is : "+this.cost_hungarian+"\n");
		this.cost_greedy = m.computeGreedyMatching(numNodes, destinationIndices);
		//s.append("The cost of matching produced by Online Greedy Algorithm is : "+this.cost_greedy+"\n");
		this.cr_grbyoff = (this.cost_greedy)/(this.cost_hungarian);
		//s.append("Competitive Ratio: greedy/hungarian : "+this.cr_grbyoff+"\n");
		
		cal = Calendar.getInstance();
		System.out.println("\nEnd: " + dateFormat.format(cal.getTime())+"\n");
		System.out.println("==========================END OF A RUN===================================\n");

		return s.toString();
	}


	public static void main(String args[]) throws NumberFormatException, IOException
	{
		//set of constant values
		BufferedWriter bw = new BufferedWriter(new FileWriter("Output.txt"));
		StringBuilder s = new StringBuilder();
		Driver driver = new Driver();
		for(int j = 0; j < args.length; j++)
		{
			int k = 1;
			while(k <= 10)
				{
					s.append(driver.generateSingleRun(Integer.parseInt(args[j]), null));
					k++;
				}
			s.append("\n");
		}
		bw.write(s.toString());
		System.out.println("Check the output files: Output.txt and Details.txt");
		bw.close();
	}
}
