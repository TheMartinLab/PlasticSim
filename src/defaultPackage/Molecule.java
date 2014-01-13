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
 * An interface that determines what a generalized molecule knows how to do:
 * Rotate
 * Translate
 * Calculate its scattering based on a given hkl
 * Output itself as a String for an Atoms .inp file
 * @author Eric D. Dill eddill@ncsu.edu
 * @author James D. Martin jdmartin@ncsu.edu
 * Copyright © 2010-2013 North Carolina State University. All rights reserved
 */
public abstract class Molecule implements Cloneable {
	protected Molecule[] surrounding;
	protected int charge;
	/**
	 * Method to rotate the molecule about an axis by an angle and an origin
	 * @param axis	The axis with which to rotate about
	 * @param origin	A quaternion representing the origin
	 * @param angle	The angle with which to rotate
	 */
	public abstract void rotate(JVector axis, JVector origin, double phi);
	
	/**
	 * Method to translate a molecule by a given amount determined by a Quaternion
	 * @param amount	The quaternion indicating the direction and magnitude with which to translate
	 */
	public abstract void translate(JVector amount);
	
	public abstract void moveTo(JVector newPos);
	
	/**
	 * toString() method to output the coordinate list of the atoms in the molecule to create an Atoms .inp file
	 * @param a	The a lattice constant
	 * @param b The b lattice constant
	 * @param c The c lattice constant
	 * @param number	The number of the atom
	 * @param aUnits	The number of units in the a direction
	 * @param bUnits	The number of units in the b direction
	 * @param cUnits	The number of units in the c direction 
	 * @return	The string to print out to the .inp file for Atoms input
	 */
	public abstract String toStringForAtoms(double a, double b, double c, int number, int aUnits, int bUnits, int cUnits);
	
	public abstract JAtom[] getLigands();
	public abstract JAtom getCenter();

	public abstract void setLigands(JAtom[] ligands);

	public abstract JVector[] getVectors();

	public abstract Object clone();

	public abstract JAtom[] getAtoms();
	public int getCharge() { return charge; }
	public void setCharge(int toSet) { charge = toSet; }
	
	public String toStringForXYZ() {
		String str = "";
		for(JAtom atom : getLigands()) {
			str += atom.toStringForXYZ();
		}
		
		return str;
	}
	
}
