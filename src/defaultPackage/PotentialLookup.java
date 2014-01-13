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

import java.util.HashMap;

public class PotentialLookup {

	protected Potential pot;
	protected double precision;
	private HashMap<Double, Double> energyLookup;
	private HashMap<Double, Double> forceLookup;
	
	public PotentialLookup(double precision, Potential pot) {
		this.pot = pot;
		this.precision = precision;
		createPotentials();
	}
	
	public boolean check(int Z1, int Z2) {
		if(Z1 == pot.Z1 && Z2 == pot.Z2)
			return true;
		if(Z1 == pot.Z2 && Z2 == pot.Z1)
			return true;
		return false;
	}
	public int getZ1() { return pot.Z1; }
	public int getZ2() { return pot.Z2; }
	public double getPotVal() { return 0; }
	public double getDistVal() { return 0; }
	public String toString() { return pot.toString(); }
	public String toStringParams() { return pot.toString(); }
	public String getType() { return "Default"; }
	
	public double lookupForce(double distance)
	{
		if(distance + precision > pot.getDistVal() * 2.5 || distance == 0)
			return 0;
		
		distance = Math.rint(distance / precision) * precision;
		
		double force=0;
		try {
			force = forceLookup.get(distance);
		} catch(NullPointerException npe) {
			force = pot.calcF(distance);
			forceLookup.put(distance, force);
		}
	
		return force;
	}
	
	public double lookupPotential(double distance)
	{
		if(distance + precision > pot.getDistVal() * 2.5 || distance == 0)
			return 0;
		
		distance = Math.rint(distance / precision) * precision;
		double potential=0;
		try {
			potential = energyLookup.get(distance);
		} catch(NullPointerException npe) {
			potential = pot.calcU(distance);
			energyLookup.put(distance, potential);
		}

		return potential;
	}
	public void createPotentials() {
		int numDataPoints = (int) (2.5 * pot.getDistVal() / precision);
		
		double maxValue = numDataPoints * precision;
		
		energyLookup = new HashMap<Double, Double>(numDataPoints);
		
		forceLookup = new HashMap<Double, Double>(numDataPoints);
		
		for(double i = 0; i <= maxValue; i += precision)
		{
			double distance = Math.rint(i / precision) * precision;
			
			double potential = pot.calcU(distance);
			
			double force = pot.calcF(distance);
			
			energyLookup.put(distance, potential);
			
			forceLookup.put(distance, force);
		}
	}
}
