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
 * @author Eric D. Dill eddill@ncsu.edu
 * @author James D. Martin jdmartin@ncsu.edu
 * Copyright © 2010-2013 North Carolina State University. All rights reserved
 */
public class JPixel implements Comparable<JPixel> {

	private JComplex[] sf;
	private JVector q;
	private double I;
	private int[] atomTypes;
	public JPixel(int numSf, JVector q, int[] atomTypes) {
		sf = new JComplex[numSf];
		this.q = (JVector) q.clone();
		setI(0);
		this.atomTypes = atomTypes;
	}
	public JPixel(int numSf) {
		sf = new JComplex[numSf];
		setI(0);
		q = new JVector();
	}
	public double getI() { return I; }
	public void setQ(JVector q) {
		this.q.i = q.i;
		this.q.j = q.j;
		this.q.k = q.k;
	}
	public void setSf(int idx, JComplex toSet) { sf[idx] = (JComplex) toSet.clone(); }
	@Override
	public int compareTo(JPixel o) {
		double len1 = q.length();
		double len2 = o.q.length();
		if(len1 < len2) { return -1; }
		else if(len1 == len2) { return 0; }
		else { return 1; }
	}
	public void setI(double i) {
		I = i;
	}
	public JComplex[] getSf() {
		return sf;
	}
	public void setSf(JComplex[] sf) {
		this.sf = sf;
	}
	public JVector getQ() {
		return q;
	}
	@Override
	public String toString() {
		String info = "Q: " + q.toString() + " I: " + I + "Scattering factors: ";
		for(int i = 0; i < atomTypes.length; i++) {
			info += atomTypes[i] + ": " + sf[i].toString();
		}
		
		return info;
	}
	public int[] getAtomTypes() {
		return atomTypes;
	}
	public void setAtomTypes(int[] atomTypes) {
		this.atomTypes = atomTypes;
	}
}
