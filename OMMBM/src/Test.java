import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;

public class Test {
	public static void main(String args[]) throws IOException
	{
		BufferedWriter bw = new BufferedWriter(new FileWriter("out"));
		SyntheticData sd = new SyntheticData();
		
	
		System.out.println("This is the beginning to your success\n");
		Element edge = new Element(10, 20, 5);
		System.out.println("The edge created is ("+edge.getX()+", "+edge.getY()+") with weight as "+edge.getWeight()+"\n");
		
		double costMatrix[][] = sd.generateSynthetic1D(5);
		System.out.println("===============================================");
		System.out.println("This is the cost matrix");
		System.out.println("===============================================");
		for(int i = 0; i < costMatrix.length; i++)
		{
			for(int j = 0; j < costMatrix[i].length; j++)
			{
				System.out.print(costMatrix[i][j]+"\t");
			}
			System.out.println("\n");
		}
		

		int V = 4;
		EdgeWeightedDigraph g = new EdgeWeightedDigraph(V);
		
		DirectedEdge e1 = new DirectedEdge(0, 1, 10);
		g.addEdge(e1);
		DirectedEdge e2 = new DirectedEdge(0, 3, 5);
		g.addEdge(e2);
		DirectedEdge e3 = new DirectedEdge(2, 1, 4);
		g.addEdge(e3);
		DirectedEdge e4 = new DirectedEdge(3, 1, 1);
		g.addEdge(e4);
		DirectedEdge e5 = new DirectedEdge(1, 2, 6);
		g.addEdge(e5);
		
	    // compute shortest paths
        DijkstraSP sp = new DijkstraSP(g, 0);
        for(int i = 0; i < g.V(); i++)
        	System.out.println("Distance to "+i+" from 0 "+sp.distTo(i)+" with path "+sp.pathTo(i));
		
		/**System.out.println("Number of vertices in the graph : "+g.V());
		System.out.println("Number of edges in the graph : "+g.E());
		System.out.println("The indegree of vertex 1 is "+g.indegree(1));
		System.out.println("The degree of vertex 0 is "+g.outdegree(0));
		
		System.out.println("\n \nThe edges in the graph are: ");
		Iterator<DirectedEdge> it = g.edges().iterator(); 
		while(it.hasNext())
		{
			DirectedEdge e = it.next();
			System.out.println(e.toString());
		}
		
		System.out.println("\n\nThe edges incident from vertex 0 are: ");
		Iterator<DirectedEdge> it1 = g.adjacentTo(0).iterator();
		while(it1.hasNext())
		{
			DirectedEdge e = it1.next();
			System.out.println(e.toString());
		}
		
		System.out.println("\n\nThe entire graph in the form of a string \n"+g.toString());
		g.reverseDirection(3, 1, 5);
		System.out.println("\n\nThe entire graph in the form of a string \n"+g.toString());
		System.out.println("The final indegree of 0, 1, 2, 3: "+g.indegree(0)+" "+g.indegree(1)+" "+g.indegree(2)+" "+g.indegree(3));
		**/
		
		/**Dijkstras dijkstra = new Dijkstras(V, 0);
		dijkstra.dsp_algorithm(g);
		**/
		bw.close();
	}
}
