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
public abstract class Potential {

	protected int Z1, Z2;
	
	public Potential() {
		this.Z1 = 0;
		this.Z2 = 0;
	}
	public Potential(int Z1, int Z2)
	{
		this.Z1 = Z1;
		this.Z2 = Z2;
	}
	
	public abstract double calcF(double r);
	public abstract double calcU(double r);
	public boolean check(int Z1, int Z2) {
		if(Z1 == this.Z1 && Z2 == this.Z2)
			return true;
		if(Z1 == this.Z2 && Z2 == this.Z1)
			return true;
		return false;
	}
	public int getZ1() { return Z1; }
	public int getZ2() { return Z2; }
	public double getPotVal() { return 0; }
	public double getDistVal() { return 0; }
	/**
	 * @return tab delimited
	 */
	public String toString() { return Z1 + "\t" + Z2; }
	/**
	 * 
	 * @return space delimited
	 */
	public String toStringParams() { return Z1 + " " + Z2; }
	public abstract String getType();
}
