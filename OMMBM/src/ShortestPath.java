/** A class to find the minimum t-net-cost augmenting path and the 
 * minimum cost
 * 
 * @author prathyush
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
	 * @return the minimum cost augmenting path
	 */
	public ArrayList<DirectedEdge> getMinCostPath()
	{
		return this.minCostPath;
	}

	public void computeMinCostPath(EdgeWeightedgraph spt, int src, int dest)
	{
		ArrayList<DirectedEdge> neighborsOfSrc = spt.adjacentTo(src);
		Iterator<DirectedEdge> it = neighborsOfSrc.iterator();
		ArrayList<DirectedEdge> temp;

		while(it.hasNext())
		{
			DirectedEdge currentEdge= it.next();
			int currentNeighbor = currentEdge.to();
			temp = new ArrayList<DirectedEdge>();
			temp.add(currentEdge);
			//System.out.println("current neighbor : "+currentNeighbor);

			if(currentNeighbor == dest)
			{
				this.minCostPath.addAll(temp);
				return;
			}

			while(spt.adjacentTo(currentNeighbor) != null)
			{
				Iterator<DirectedEdge> it1 = spt.adjacentTo(currentNeighbor).iterator();
				currentEdge = it1.next();
				currentNeighbor = currentEdge.to();
				temp.add(currentEdge);
				//System.out.println("current neighbor : "+currentNeighbor);
				if(currentNeighbor == dest)
				{
					this.minCostPath.addAll(temp);
					return;
				}
			}
		}
	}



}
