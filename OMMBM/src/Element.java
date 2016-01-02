/* Element represents an Edge in the graph.
 */
public class Element {

	private int x;
	private int y;
	private double weight;

	/**
	 * Creates a new element
	 * 
	 * @param x
	 *            Index of edge in Set A
	 * @param y
	 *            Index of edge in Set B
	 * @param weight
	 *            Weight/Cost of edge
	 */
	public Element(int x, int y, double weight) {
		this.setX(x);
		this.setY(y);
		this.setWeight(weight);
	}

	/**
	 * Getters and Setters
	 */

	public double getWeight() {
		return weight;
	}

	public void setWeight(double weight) {
		this.weight = weight;
	}

	public int getY() {
		return y;
	}

	public void setY(int y) {
		this.y = y;
	}

	public int getX() {
		return x;
	}

	public void setX(int x) {
		this.x = x;
	}

	/**
	 * Returns a string representation of an Element
	 */
	public String toString() {
		return this.getWeight() + ": (" + this.getX() + ", " + this.getY()
				+ ")";
	}

}
