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
import java.io.Serializable;

public class HookePotential extends Potential implements Cloneable, Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 7932373660939712090L;
	/** Spring constant **/
	protected double k;
	/** Equilibrium distance **/
	protected double r0;
	
	public HookePotential(int Z1, int Z2, double r0, double k)
	{
		super(Z1, Z2);
		this.k = k;
		this.r0 = r0;
	}
	@Override
	public double calcF(double r) { return - k *(r-r0); }
	@Override
	public double calcU(double r) { return -.5*Math.pow(r-r0, 2) * k; }
	@Override
	public Object clone() { return new HookePotential(Z1, Z2, k, r0); }
	@Override
	public double getPotVal() { return k; }
	@Override
	public double getDistVal() { return r0; }
	@Override	
	public String toString() { return getType() + "\t" + super.toString() + "\t" + r0 + "\t" + k + "\t"; }
	@Override
	public String toStringParams() { return getType() + " " + super.toStringParams() + " " + r0 + " " + k; }
	@Override
	public String getType() { return "Hooke"; }
}
