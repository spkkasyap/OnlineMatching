import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;


/** Computes the offline matching and online matching using 
   Dijkstra's algorithm and dual weights. Calculates the Competitive Ratio
   to evaluate the performance of our online algorithm.

    @author prathyush
 **/

public class Matching {
	private double[][] costMatrix;
	
	/** A function to compute the cost matrix from the datasource
	 * 
	 */
	public void generateCostMatrix(int numSetA)
	{
		SyntheticData sd = new SyntheticData();
		this.costMatrix = sd.generateSynthetic1D(numSetA);
	}
	
	/** A function to print the cost matrix 
	 * 
	 */
	public void printCostMatrix()
	{
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
	}

	/**
	 * Generates a randomized list of request or destination indices for the
	 * online algorithm to utilize.
	 * 
	 * @param numSetA
	 *            The number of taxis (nodes in set A)
	 * @return A randomized ArrayList of destination indices: numSetA <= i <
	 *         numSetA*2
	 */
	public ArrayList<Integer> permuteDestinations(int numSetA) {

		ArrayList<Integer> destinationIndices = new ArrayList<Integer>();

		for (int i = numSetA; i < numSetA * 2; i++) {
			destinationIndices.add(i);
		}

		Collections.shuffle(destinationIndices);

		return destinationIndices;
	}



	/**
	 * Helper method that checks if two DirectedEdges are equal to each other.
	 * To, from and weight attributes must equal in both objects.
	 * 
	 * @param edge1
	 *            First directed edge
	 * @param edge2
	 *            Second directed edge
	 * @return True if all attributes are equal in both edges. False if either
	 *         edge is null or attributes don't match.
	 */
	public boolean edgeEquals(DirectedEdge edge1, DirectedEdge edge2) {

		if (edge1.equals(null) || edge2.equals(null)) {
			return false;
		}

		if (edge1.to() == edge2.to() && edge1.from() == edge2.from()
				&& edge1.weight() == edge2.weight()) {
			return true;
		}

		return false;
	}


	/**
	 * Helper method that checks if an ArrayList of DirectedEdges contains an
	 * edge.
	 * 
	 * @param list
	 *            The list of DirectedEdges that may or may not contain the edge
	 * @param edge
	 *            The edge to check if it is contained in the list
	 * @return True if the list contains edge, false otherwise.
	 */
	public boolean containsEdge(ArrayList<DirectedEdge> list, DirectedEdge edge) {
		for (DirectedEdge e : list) {
			if (edgeEquals(e, edge)) {
				return true;
			}
		}

		return false;
	}


	/**
	 * Calculates the total net cost of the matching computed by the selected
	 * algorithm
	 * 
	 * @param matching
	 *            The edges in the matching from running the selected algorithm
	 * @return Total net cost of the matching
	 */
	public double calculateTotalCost(ArrayList<DirectedEdge> matching) {
		double testCost = 0;
		for (DirectedEdge test : matching) {
			testCost += Math.abs(test.weight());
		}
		return testCost;
	}

	/**
	 * Builds the Element ArrayList when given a cost matrix
	 * 
	 * @param costMatrix
	 *            A simple 2-D array representing the costs between edges
	 * @return An ArrayList of Elements used to create DirectedEdges
	 */
	public static ArrayList<Element> buildElementMatrix(double[][] costMatrix) {

		ArrayList<Element> elementList = new ArrayList<Element>();

		for (int i = 0; i < costMatrix.length; i++) {
			for (int j = 0; j < costMatrix[0].length; j++) {
				Element element = new Element(i, j + costMatrix.length, costMatrix[i][j]);
				elementList.add(element);
			}
		}
		return elementList;
	}

	/**
	 * Builds the EdgeWeightedDigraph when given a cost matrix
	 * 
	 * @param costMatrix
	 *            A simple 2-D array representing the costs between edges
	 * @return The EdgeWeightedDigraph; simply a bipartite graph with all
	 *         vertices, edges and weights associated in one data structure
	 */
	public static EdgeWeightedgraph constructDigraphFromMatrix(double[][] costMatrix) {
		EdgeWeightedgraph diGraph = new EdgeWeightedgraph(costMatrix.length
				+ costMatrix[0].length);
		
		ArrayList<Element> elements = buildElementMatrix(costMatrix);

		for (Element ele : elements) {
			DirectedEdge edge = null;
			if (ele.getWeight() < 0) {
				edge = new DirectedEdge(ele.getY(), ele.getX(), ele.getWeight());
			} else {
				edge = new DirectedEdge(ele.getX(), ele.getY(), ele.getWeight());
			}
			diGraph.addEdge(edge);
		}
		return diGraph;
	}
	

	/**
	 * Computes the smallest cost matching using the Bellman ford algorithm in
	 * the offline setting.
	 * 
	 * @param filename
	 *            Name of dataset file.
	 * @param numSetA
	 *            The number of taxis (nodes in set A)
	 * 
	 * @return The total net cost of the best offline matching
	 */
	public double computeOfflineMatching(int numSetA) {
		// Final offline matching ArrayList of directed edges
		ArrayList<DirectedEdge> offlineMatching = new ArrayList<DirectedEdge>();

		// Temporary matching ArrayList of Set B values for internal management
		ArrayList<Integer> matching = new ArrayList<Integer>();

		// ArrayList that stores all negative cycle indices for repeated
		// processing
		ArrayList<Integer> negativeCycleIndex = new ArrayList<Integer>();

		// Index of the current source node being processed
		int index = 0;

		// Boolean value that is set to true
		boolean runNegativeCycleIndices = false;

		// Construct a DiGraph from the original costmatrix
		EdgeWeightedgraph original = constructDigraphFromMatrix(costMatrix);

		ArrayList<Integer> sourceIndices = new ArrayList<Integer>();
		for (int i = 0; i < numSetA; i++) {
			sourceIndices.add(i);
		}

		// Randomizes the source node indices chosen for augmentation
		Collections.shuffle(sourceIndices);

		/*
		 * Core of the algorithm
		 * 
		 * Algorithm repeats until matching is perfect. Source node index is
		 * picked and BellmanFord is ran on that node. This computes all the
		 * paths and their costs from the source node to every target node. All
		 * paths are processed to choose the minimum cost and to also check for
		 * negative cycles. If a negative cycle occurs on the node, then add it
		 * to the negativeCycleIndex ArrayList and move on to the next source
		 * node. Else, the minimum cost path is chosen. The best path
		 * (augmenting or direct path) is added to an ArrayList bestPath and a
		 * new DiGraph is constructed with the best cost matching and the
		 * previous matchings. This new edge matching is added to the matchings
		 * ArrayList and the original DiGraph is updated with this new DiGraph.
		 * The process is repeated.
		 * 
		 * If negative cycles occurred during the processing, an additional
		 * iteration over the negativeCycleIndex ArrayList is ran to process
		 * these indices and find their proper matchings.
		 */
		while (matching.size() < numSetA) {

			if (index == numSetA) { // Last source index to process?
				if (negativeCycleIndex.size() > 0) { // Negative cycles to
														// process?
					index = negativeCycleIndex.remove(0);
					runNegativeCycleIndices = true;
				} else {
					break; // All processing complete
				}
			}

			int source = sourceIndices.get(index);

			// Run BellmanFord algorithm on source index
			BellmanFordSP sp = new BellmanFordSP(original, source);

			Iterator<DirectedEdge> iter = null;
			ArrayList<DirectedEdge> bestPath = new ArrayList<DirectedEdge>();

			double minPath = Double.MAX_VALUE;

			// Obtain minimum cost path
			for (int v = numSetA; v < original.V(); v++) {

				// If vertex is already in the matching, skip it
				if (matching.contains(v)) {
					continue;
				}

				// Check if vertex causes a negative cycle
				if (sp.hasNegativeCycle()) {
					negativeCycleIndex.add(index);
					index++;
					break;
				}

				// Check if a path exists from source vertex to destination
				// vertex v
				if (sp.hasPathTo(v)) {
					if (sp.distTo(v) < minPath) {
						minPath = sp.distTo(v);
						iter = sp.pathTo(v).iterator();
					}
				}
			}

			// Negative cycle detected, don't process current index
			if (iter == null) {
				if (runNegativeCycleIndices) {
					index = numSetA;
				}
				continue;
			}

			while (iter.hasNext()) {
				bestPath.add(iter.next());
			}

			// New iteration DiGraph that will have updated edges
			EdgeWeightedgraph nextIterationGraph = new EdgeWeightedgraph(numSetA
					+ costMatrix[0].length);

			// Update matchings and new DiGraph
			for (DirectedEdge e : original.edges()) {
				if (bestPath.contains(e)) {
					DirectedEdge edgeToAdd = new DirectedEdge(e.to(), e.from(), -1.0 * e.weight());
					nextIterationGraph.addEdge(edgeToAdd);
					if (edgeToAdd.weight() <= 0) {
						if (!matching.contains(e.to())) {
							matching.add(e.to());
						}
					}
				} else {
					nextIterationGraph.addEdge(e);
				}
			}

			original = nextIterationGraph;

			// Increment index if not processing negative cycles, else set index
			// to last source node
			if (!runNegativeCycleIndices) {
				index++;
			} else {
				index = numSetA;
			}
		}

		// Add all matched edges to offline matching with proper weight
		for (DirectedEdge edge : original.edges()) {
			if (edge.weight() < 0) {
				offlineMatching.add(new DirectedEdge(edge.to(), edge.from(), -1.0 * edge.weight()));
			}
		}

		double totalCost = calculateTotalCost(offlineMatching);

		return totalCost;
	}

	
	
	/**
	 * Computes the smallest cost matching using the Hungarian algorithm. For
	 * verification purposes.
	 * 
	 * @return The net total cost of the optimal matching found by executing
	 *         Hungarian algorithm.
	 */
	public double verifyHungarian() {
		HungarianAlgorithm test = new HungarianAlgorithm(costMatrix);

		int[] tester = test.execute();
		double totalCost = 0;
		for (int i = 0; i < tester.length; i++) {
			totalCost += costMatrix[i][tester[i]];
		}
		return totalCost;
	}
	
}
