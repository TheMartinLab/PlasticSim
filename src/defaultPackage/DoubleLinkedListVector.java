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
 * A doubly linked list class.  This class can insert and remove from 
 * the end of the list.  This class uses the TwoNode class.
 * @author Eric D. Dill eddill@ncsu.edu
 * @author James D. Martin jdmartin@ncsu.edu
 * Copyright © 2010-2013 North Carolina State University. All rights reserved
 */
public class DoubleLinkedListVector {

	/** Pointers to the head and tail of the list */
	private TwoNodeVector first;
	private TwoNodeVector last;
	
	protected int length;
	
	private final int maxLength;
	
	/**
	 * Constructor to initialize the head and tail of the list to each other with null data in each.
	 */
	public DoubleLinkedListVector(int maxLength)
	{
		first = new TwoNodeVector(0, 0, null, null);
		last = new TwoNodeVector(0, 0, null, null);
		first.setNext(last);
		last.setPrev(first);
		length = 0;
		this.maxLength = maxLength;
	}
	public void clear() {
		while(first.next!= last)
			removeHead();
	}
	public boolean hasNext()
	{
		if(first.getNext() != last)
			return true;
		
		return false;
	}
	public void insert(TwoNodeVector aNode)
	{
		TwoNodeVector temp = findWhereToInsert(aNode.value);
		
		if(length == maxLength && temp.getNext() != last)
		{
			insertBefore(temp, aNode);
			removeLast();
		}
		
		if(length < maxLength)
		{
			insertBefore(temp, aNode);
		}
	}
	
	
	public TwoNodeVector findWhereToInsert(double angle)
	{
		TwoNodeVector temp = first.getNext();
		
		while(temp != last && angle > temp.value)
		{
			temp = temp.getNext();
		}
		
		return temp;
	}
	/**
	 * Method to insert an existing node at the head of the list
	 * @param toInsert	the existing node to insert at the head of the list.
	 */
	public void insertHead(TwoNodeVector toInsert)
	{
		toInsert.setPrev(first);
		toInsert.setNext(first.getNext());
		first.getNext().setPrev(toInsert);
		first.setNext(toInsert);
		length++;
	}
	/**
	 * Method to insert a Node at the tail of the list
	 * @param data1 first object to be contained in the node.
	 * @param data2 second object to be contained in the node.
	 */
	public void insertTail(double value)
	{
		last.getPrev().setNext(new TwoNodeVector(value, 0, last.getPrev(), last));
		last.setPrev(last.getPrev().getNext());
		length++;
	}
	/**
	 * Method to insert a node after a specified node in the list.
	 * @param pos		Node to insert after.
	 * @param toInsert	Node to insert.
	 */
	public void insertAfter(TwoNodeVector pos, TwoNodeVector toInsert)
	{
		toInsert.setPrev(pos);
		toInsert.setNext(pos.getNext());
		pos.getNext().setPrev(toInsert);
		pos.setNext(toInsert);
		length++;
	}
	/**
	 * Method to insert a node before the specified node.
	 * @param pos		Node to insert before.
	 * @param toInsert	Node to insert.
	 */
	public void insertBefore(TwoNodeVector pos, TwoNodeVector toInsert)
	{
		insertAfter(pos.getPrev(), toInsert);
	}
	/**
	 * Method to remove a specified node from the list.
	 * @param toRemove	Node to remove.
	 * @return			Node that was removed.
	 */
	public TwoNodeVector remove(TwoNodeVector toRemove)
	{
		toRemove.getPrev().setNext(toRemove.getNext());
		toRemove.getNext().setPrev(toRemove.getPrev());
		length--;
		return toRemove;
	}
	
	public TwoNodeVector removeHead()
	{
		TwoNodeVector temp = first.getNext();
		
		temp.getNext().setPrev(first);
		
		first.setNext(temp.getNext());
		
		temp.setNext(null);
		temp.setPrev(null);
		
		length--;
		
		return temp;
	}
	public TwoNodeVector removeLast()
	{
		TwoNodeVector temp = last.getPrev();
		last.setPrev(last.getPrev().getPrev());
		last.getPrev().setNext(last);
		
		temp.setNext(null);
		temp.setPrev(null);
		
		length--;
		
		return temp;
	}
	/**
	 * Method to return the list as a comma-delimited string.
	 * @return	the list as a string.
	 */
	public String toString()
	{
		String s = "";
		TwoNodeVector temp = first.getNext();
		while(temp != last)
		{
			s += temp.toString();
			//s += temp.getData() + " " + temp.getCount() + "\n";
			temp = temp.getNext();
		}
		return s;
	}
	
	public TwoNodeVector getLast() { return last.getPrev(); }
	
	public boolean checkLength() { return length == maxLength; }
}
