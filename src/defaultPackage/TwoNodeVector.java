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
public class TwoNodeVector{
	
	public TwoNodeVector prev;	// pointer to the previous node.
	public TwoNodeVector next;	// pointer to the next node.
	protected double value;	// String contained in the node.
	protected int key;
	private int key2;
	
	/**
	 * Constructor to initialize a node.
	 * @param d			The value to be stored in the list.
	 * @param key		The sorting value.
	 * @param previous	reference to the previous node.
	 * @param next		reference to the next node.
	 */
	public TwoNodeVector(double value, int key, TwoNodeVector previous, TwoNodeVector nextNode)
	{
		next = nextNode;
		prev = previous;
		this.value = value;
		this.key = key;
	}
	/**
	 * Setter method for the reference to the previous node.
	 * @param previous	reference to the previous node.
	 */
	public void setPrev(TwoNodeVector previous) { prev = previous; }
	/**
	 * Setter method for the reference to the next node.
	 * @param nextNode reference to the next node.
	 */
	public void setNext(TwoNodeVector nextNode) { next = nextNode; }
	/**
	 * Getter method for the previous node.
	 * @return the previous node
	 */
	public TwoNodeVector getPrev() {  return prev; }
	/**
	 * Getter method for the next node.
	 * @return the next node.
	 */
	public TwoNodeVector getNext() { return next; }
	/**
	 * Getter method for the Object in the node.
	 * @return the string in the node.
	 */
	public double getValue() { return value; }
	
	public int getKey() { return key; }

	/**
	 * Return the String and counter variable of the node as a 
	 * comma-delimited string.
	 */
	public String toString() { return key + "\t" + value; }
	public int getKey2() {
		return key2;
	}
	public void setKey2(int key2) {
		this.key2 = key2;
	}
}
