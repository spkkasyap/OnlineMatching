/** A class to find the minimum t-net-cost augmenting path and the 
 * minimum cost
 */

import java.util.ArrayList;
import java.util.Iterator;

public class ShortestPath {
	private double minCost;
	private ArrayList<DirectedEdge> minCostPath;

	/** The constructor to initialize the private variables of the class
	 */
	public ShortestPath()
	{
		minCost = Double.MAX_VALUE;
		minCostPath = new ArrayList<DirectedEdge>();
	}

	/** A function to set the minimum cost to a particular value
	 * 
	 * @param x
	 */
	public void setMinCost(int x)
	{
		this.minCost = x;
	}

	/** A function to add an edge to the minimum cost 
	 * augmenting path 
	 * 
	 * @param e
	 */
	public void addToMinCostPath(DirectedEdge e)
	{
		this.minCostPath.add(e);
	}

	/** A function to set the Shortest path to a particular path
	 * 
	 * @param path
	 */
	public void setMinCostPath(ArrayList<DirectedEdge> path)
	{
		Iterator<DirectedEdge> it = path.iterator();
		while(it.hasNext())
		{
			this.minCostPath.add(it.next());
		}
	}
	
	/** A function to get the minimum cost 
	 * 
	 * @return minimum cost of the shortest path
	 */
	public double getMinCost()
	{
		return this.minCost;
	}
	
	/** A function to get the shortest path
	 * 
	 * @return the minimun 
	 */
	public ArrayList<DirectedEdge> getMinCostPath()
	{
		return this.minCostPath;
	}


}
