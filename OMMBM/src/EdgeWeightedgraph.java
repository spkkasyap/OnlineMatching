import java.util.ArrayList;
import java.util.Iterator;

/** 
 This is a class that holds an edge weighted directed graph
 @author prathyush
 **/
public class EdgeWeightedgraph {

	private int V; 	//number of vertices in the graph
	private int E;	//number of edges in the graph
	private ArrayList<DirectedEdge>[] adj;	//adjacency list of the vertices
	private int [] indegree; //in degree of the vertices of the graph

	/** A constructor to create an empty directed graph with V vertices
	 *
	 * @param V (integer)
	 */


	public EdgeWeightedgraph(int V)
	{
		this.V = V;
		this.E = 0;

		this.indegree = new int[V];
		adj = (ArrayList<DirectedEdge>[]) new ArrayList[V];
		for (int i = 0; i < V; i++)
			adj[i] = null;
	}

	/** A function to return the number of vertices in the graph
	 * 
	 * @return V; number of vertices in the graph
	 */

	public int V()
	{
		return V;
	}

	/** A function to return the number of edges in the graph
	 * 
	 * @return E; number of edges in the graph
	 */
	public int E()
	{
		return E;
	}


	/** A function to add a directed edge to the directed graph 
	 * 
	 * @param edge
	 */
	public void addEdge(DirectedEdge edge)
	{
		int v = edge.from();
		int w = edge.to();
		if(adj[v] == null)
			adj[v] = new ArrayList<DirectedEdge>();
		adj[v].add(edge);
		indegree[w]++;
		E++;
	}

	/** A function to return the indegree of the vertex given as the parameter
	 * 
	 * @param v
	 * @return
	 */
	public int indegree(int v)
	{
		return indegree[v];
	}

	/** A function to return the outdegree of the vertex given as the parameter
	 * 
	 * @param v
	 * @return
	 */
	public int outdegree(int v)
	{
		return adj[v].size();
	}

	/** A function to return all the edges incident from given vertex v
	 * 
	 * @param v
	 * @return adjacency list of v.
	 */
	public ArrayList<DirectedEdge> adjacentTo(int v)
	{
		return adj[v];
	}

	/** A function to return all the directed edges in the graph
	 * 
	 * @return allEdges ; all directed edges in the graph
	 */
	public ArrayList<DirectedEdge> edges()
	{
		ArrayList<DirectedEdge> allEdges = new ArrayList<DirectedEdge>();
		for(int i = 0; i < V; i++)
		{
			if(adj[i] != null)
				allEdges.addAll(adj[i]);
		}
		return allEdges;
	}

	/** A function to remove edges incident on a particular vertex
	 * 
	 */
	public void removeInto(int u)
	{
		ArrayList<DirectedEdge> toRemove = new ArrayList<DirectedEdge>();
		for(int i =0; i < this.V(); i++)
		{
			if(this.adjacentTo(i) != null)
			{
				Iterator<DirectedEdge> it = this.adjacentTo(i).iterator();
				while(it.hasNext())
				{
					DirectedEdge e = it.next();
					if(e.to() == u)
					{
						toRemove.add(e);
					}
				}
			}
		}
		Iterator<DirectedEdge> it1 = toRemove.iterator();
		while(it1.hasNext())
		{
			DirectedEdge e1 = it1.next();
			this.removeEdge(e1);
		}
			
	}

	/** A function to remove an edge from the graph. This method requires to 
	 * keep track of the effects of removing an edge on other parameters like 
	 * E, indegree[].
	 * 
	 * @param e
	 */
	public void removeEdge(DirectedEdge e)
	{
		int v = e.from();
		int w = e.to();
		this.adj[v].remove(e);
		this.indegree[w]--;
		this.E--;
	}

	/** A function to reverse an edge direction as well to change its weight. Given the initial direction
	 * to be from v to w, we first remove the edge v->w and then, add the edge 
	 * w->v with new weight 'weight'. This function uses 'removeEdge()' as a
	 * sub module.
	 * @param v
	 * @param w
	 * @param weight
	 */
	public void reverseDirection(int v, int w, double weight)
	{
		DirectedEdge e = null;
		Iterator<DirectedEdge> it = this.adjacentTo(v).iterator();
		while(it.hasNext())
		{
			e = it.next();
			int f = e.to(); 
			System.out.println(f);
			if(f == w)
			{
				break;
			}
		}
		this.removeEdge(e);
		DirectedEdge e1 = new DirectedEdge(w, v, weight);
		this.addEdge(e1);
		//System.out.println("E at this point "+E);
	}

	/** A function to write down the entire directed graph as a string
	 * 
	 * @return string form of the edge weighted directed graph
	 */
	public String toString()
	{
		StringBuilder s = new StringBuilder();
		s.append("V: "+V+" E :"+E+"\n");
		Iterator<DirectedEdge> it = this.edges().iterator();
		while(it.hasNext())
		{
			DirectedEdge e = it.next();
			s.append(e.toString()+"\n");
		}
		return s.toString();
	}

}
