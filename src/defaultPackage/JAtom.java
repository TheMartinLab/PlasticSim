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
import java.io.Serializable;
/**
 * @author Eric D. Dill eddill@ncsu.edu
 * @author James D. Martin jdmartin@ncsu.edu
 * Copyright © 2010-2013 North Carolina State University. All rights reserved
 */
import java.util.Observable;

/**
 * This class represents an JAtom.  An JAtom knows its position in (x, y, z) coordinates and its Z.
 * @author Eric Dill
 * 			eddill@ncsu.edu
 *
 */
public class JAtom implements Cloneable, Serializable
{

	/** The abbreviation of the JAtom: Carbon -> C */
	
	private static final long serialVersionUID = 8225076374777785538L;

	/** The Z value of the JAtom: Carbon -> 6 */
	protected int Z;
	
	/** The position in 3 space of the JAtom */
	private JVector position;
	
	/** The total number of atoms in the system */
	protected static int numberOfAtoms;
	
	/** The identifier of the JAtom object */
	protected int atomID;

	/** 
	 * Constructor, initialize the atom based on a vector and a name
	 * @param Z		Number of electrons in the atom
	 * @param v1	The vector containing the atomic coordinates
	 */
	public JAtom(int Z, JVector v1)
	{
		setPosition((JVector) v1.clone());
		
		this.Z = Z;
		
		atomID = numberOfAtoms++;
	}
	
	/**
	 * Shift the position of the atom by the JVector amount
	 * @param amount	Amount to shift the atomic center by
	 */
	protected void translate(JVector amount)
	{
		getPosition().add_inPlace(amount);
	}
	
	/**
	 * Getter method for number of electrons in the atom
	 * @return	The atomic abbreviation, eg "H" or "He"
	 */
	public int getZ() { return Z; }

	public JVector getPosition() { return position; }
	
	public Object clone()
	{
		return new JAtom(Z, (JVector) getPosition().clone()); 
	}

	/**
	 * toString() override for the JAtom class.  Output: "Z + (i, j, k)"
	 */
	public String toString() { return Z + "," + getPosition(); }

	/**
	 * toString() method to output the coordinate list of the atom to create an Atoms .inp file
	 * @param a	The a lattice constant
	 * @param b The b lattice constant
	 * @param c The c lattice constant
	 * @param number	The number of the atom
	 * @param aUnits	The number of units in the a direction
	 * @param bUnits	The number of units in the b direction
	 * @param cUnits	The number of units in the c direction 
	 * @return	The string to print out to the .inp file for Atoms input
	 */
	public String toStringForAtoms(double a, double b, double c, int number, double aUnits, double bUnits, double cUnits)
	{
		double i = getPosition().i;
		
		double j = getPosition().j;
		
		double k = getPosition().k;
		
		return Z + "\t" + i/(a*aUnits) + "\t" + j/(b*bUnits) + "\t" + k/(c*cUnits) + "\t" + Z + "\t" + number; 
	}

	public String toStringForXYZ()
	{
		double i = getPosition().i;
		
		double j = getPosition().j;
		
		double k = getPosition().k;
		
		return Z + "\t" + i + "\t" + j + "\t" + k + "\n";
	}
	/**
	 * Multiply the vectorial position of the atom by the passed scalar amount
	 * @param c	The scalar amount to multiply the atom by
	 */
	public void multiply(double c) { setPosition(JVector.multiply(getPosition(), c)); }
	public void setNewPos(JVector newPos) { setPosition(newPos); }
	public int getAtomID() { return atomID; }
	public void setAtomID(int atomID) { this.atomID = atomID; }

	public void setPosition(JVector position) {
		this.position = position;
	}
}
