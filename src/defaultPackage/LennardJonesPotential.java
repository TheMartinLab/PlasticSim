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

public class LennardJonesPotential extends Potential implements Cloneable, Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 8758121726760837622L;
	/* Instance variables */
	/** The depth of the potential well at the equilibrim distance**/
	protected double epsilon;
	/** The equilibrium distance **/
	protected double sigma, rMin;
	
	
	public LennardJonesPotential() {
		super();
		epsilon = 0;
		sigma = 0;
	}
	/**
	 * Constructor.
	 * @param sigma	The ideal distance between the two atoms, this corresponds to a minimum in the energy well
	 * @param epsilon	The depth of the potential well
	 */
	public LennardJonesPotential(int Z1, int Z2, double sigma, double epsilon)
	{
		super(Z1, Z2);
		this.sigma = sigma;
		rMin = sigma;
		this.epsilon = epsilon;
	}
	@Override
	/**
	 * Method to calculate the energy of interaction between two atoms separated by a distance r
	 * @param r	The distance between the two atoms
	 * @return	The force
	 */
	public double calcF(double r) 
	{
		//if(r < rMin) { return -((24*epsilon*Math.pow(sigma, 6)/Math.pow(r, 7)) / (2*Math.pow(sigma, 6)/Math.pow(r, 6) - 1)); }
		return (24*epsilon*Math.pow(sigma, 6)/Math.pow(r, 7)) * (2*Math.pow(sigma, 6)/Math.pow(r, 6) - 1);
		//return (24*epsilon*Math.pow(sigma, 6)/Math.pow(r, 7)) / (2*Math.pow(sigma, 6)/Math.pow(r, 6) - 1);
		//return (48*epsilon*Math.pow(sigma, 12)/Math.pow(r, 13));
	}
	
	@Override
	/**
	 * Method to calculate the energy of interaction between two atoms separated by a distance r
	 * @param r	The distance between the two atoms
	 * @return	The potential
	 */
	public double calcU(double r) 
	{
		//if(r < sigma) { return -((4*epsilon*Math.pow(sigma, 6)/Math.pow(r, 6)) / (Math.pow(sigma, 6)/Math.pow(r, 6) - 1)); }
		return (4*epsilon*Math.pow(sigma, 6)/Math.pow(r, 6)) * (Math.pow(sigma, 6)/Math.pow(r, 6) - 1);
		//return (4*epsilon*Math.pow(sigma, 6)/Math.pow(r, 6)) / (Math.pow(sigma, 6)/Math.pow(r, 6) - 1);
		//return (4*epsilon*Math.pow(sigma, 12)/Math.pow(r, 12));
	}
	@Override
	public Object clone() { return new LennardJonesPotential(Z1, Z2, sigma, epsilon); }
	@Override
	public double getPotVal() { return epsilon; }
	@Override
	public double getDistVal() { return sigma; }
	@Override
	public String toString() { return getType() + "\t" + super.toString() + "\t" + sigma + "\t" + epsilon + "\t"; }
	@Override
	public String toStringParams() { return getType() + " " + super.toStringParams() + " " + sigma + " " + epsilon; }
	@Override
	public String getType() { return "Lennard-Jones (6-12)"; }
	public void test() {
		for(double d = 3; d < 7; d+=.01) {
			System.out.println(d + "\t" + calcF(d));
		}
	}
	public static void main(String[] args) {
		LennardJonesPotential ljp = new LennardJonesPotential(35, 35, 3.385415128933289, .12);
		ljp.test();
	}
}
