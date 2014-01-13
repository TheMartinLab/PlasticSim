/*******************************************************************************
 * Copyright (c) 2013 Eric Dill -- eddill@ncsu.edu. North Carolina State University. All rights reserved.
 * This program and the accompanying materials
 * are made available under the terms of the GNU Public License v3.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors:
 *     Eric Dill -- eddill@ncsu.edu - initial API and implementation
 ******************************************************************************/
package defaultPackage;
/**
 * A Node class that contains a String and an integer telling how many times the
 * String occurs and a reference to both the node before and the node after. 
 * This class is meant to be used with the DoubleLinkedList class.
 * @author Eric D. Dill eddill@ncsu.edu
 * @author James D. Martin jdmartin@ncsu.edu
 * Copyright © 2010-2013 North Carolina State University. All rights reserved
 */
public class TwoNodeCBr4{
	
	private TwoNodeCBr4 prev;	// pointer to the previous node.
	private TwoNodeCBr4 next;	// pointer to the next node.
	protected IdealTetrahedron value;	// String contained in the node.
	protected double key;
	
	/**
	 * Constructor to initialize a node.
	 * @param d			The value to be stored in the list.
	 * @param key		The sorting value.
	 * @param previous	reference to the previous node.
	 * @param next		reference to the next node.
	 */
	public TwoNodeCBr4(IdealTetrahedron value, double key, TwoNodeCBr4 previous, TwoNodeCBr4 nextNode)
	{
		next = nextNode;
		prev = previous;
		this.value = value;
		this.key = key;
	}
	
	public TwoNodeCBr4(IdealTetrahedron value, double key) {
		this(value, key, null, null);
	}
	/**
	 * Setter method for the reference to the previous node.
	 * @param previous	reference to the previous node.
	 */
	public void setPrev(TwoNodeCBr4 previous) { prev = previous; }
	/**
	 * Setter method for the reference to the next node.
	 * @param nextNode reference to the next node.
	 */
	public void setNext(TwoNodeCBr4 nextNode) { next = nextNode; }
	/**
	 * Getter method for the previous node.
	 * @return the previous node
	 */
	public TwoNodeCBr4 getPrev() {  return prev; }
	/**
	 * Getter method for the next node.
	 * @return the next node.
	 */
	public TwoNodeCBr4 getNext() { return next; }
	/**
	 * Getter method for the Object in the node.
	 * @return the string in the node.
	 */
	public IdealTetrahedron getValue() { return value; }
	
	public double getKey() { return key; }

	/**
	 * Return the String and counter variable of the node as a 
	 * comma-delimited string.
	 */
	public String toString() { return key + "\t" + value.toString(); }
}
