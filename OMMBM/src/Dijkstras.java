import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

/** This is a class designed to incorporate the algorithm 
 * and some of its outputs for our use in online matching
 * 
 * @author prathyush
 *
 */
public class Dijkstras {
	private double[] distance;
	private boolean[] visited;
	private HashSet<Integer> unvisited;
	private final int source;

	/** Constructor for the Dijkstra's algorithm which initializes all the required 
	 * values of various parameters
	 * @param V
	 * @param src
	 */
	public Dijkstras(int V, int src)
	{	
		this.source = src;
		this.visited = new boolean[V];
		this.unvisited = new HashSet<Integer>();
		this.distance = new double[V];

		for(int i = 0; i < V; i++) //initialization  
		{
			this.visited[i] = false; //set the visited[.] for every vertex . as false initially
			unvisited.add(i);
			if(i != source)
				this.distance[i] = Double.MAX_VALUE; //set the distance of every vertex to source as the maximum value initially
			else
				this.distance[i] = 0.0; //distance of source to itself is set as 0
		}
	}

	/** A function to find out the vertex which is closest to the source vertex
	 * from the unvisited vertices and return it
	 * @return min_distant; the minimum distant vertex from source
	 */
	public int minDistant()
	{
		int min_distant = 0;
		double min_distance =  Double.MAX_VALUE;

		Iterator<Integer> it = this.unvisited.iterator();
		while(it.hasNext())
		{
			int v = it.next();
			if(this.distance[v] < min_distance)
			{
				min_distant = v;
				min_distance = this.distance[v];
			}
		}
		return min_distant;
	}

	/** A function to display the final shortest distances from a given
	 * source vertex to every other vertex in the graph.
	 */
	public void printSolution()
	{
		System.out.println("\n===================\nSOURCE: "+this.source+"\n==================="
				+ "\nVERTEX \t DISTANCE FROM SOURCE\n===================================");
		for(int i = 0; i < this.distance.length; i++)
		{
			System.out.println(i+"\t"+this.distance[i]);
		}
	}

	/** The main algorithm which finds the shortest path from a single source vertex 
	 * to every other vertex. The shortest path to each vertex from the source is recorded
	 * in the distance[.] for every vertex .
	 * @param g ; the edge weighted directed graph
	 */
	public void dsp_algorithm(EdgeWeightedgraph g)
	{
		//double round_min = Double.MAX_VALUE;
		//double overall_min = Double.MAX_VALUE;
		//Graph shortestPathTree  = new Graph();
		EdgeWeightedgraph spt = new EdgeWeightedgraph(g.V());
		
		for(int count = 0; count < g.V(); count++)
		{
			int u = this.minDistant();
			this.unvisited.remove(u);
			this.visited[u] = true;
			
			Iterator<DirectedEdge> it = g.adjacentTo(u).iterator();
			while(it.hasNext())
			{
				DirectedEdge e = it.next();
				int v = e.to();
				double weight = e.weight();

				if(this.visited[v] == false && this.distance[v] > weight+this.distance[u])
				{
					DirectedEdge e1 = new DirectedEdge(u, v, weight);
					if(spt.indegree(v) > 0)
						spt.removeInto(v);
					spt.addEdge(e1);
					this.distance[v] = weight+this.distance[u];
				}
			}
			
		}
		
		System.out.println("\n\nShortest Path Tree :"+spt.toString());
		this.printSolution();
		
		for(int i = 0; i < this.visited.length; i++)
			System.out.println(visited[i]);
	}
}
