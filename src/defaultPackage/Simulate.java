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
/**
 * @author Eric D. Dill eddill@ncsu.edu
 * @author James D. Martin jdmartin@ncsu.edu
 * Copyright © 2010-2013 North Carolina State University. All rights reserved
 */
package defaultPackage;
import io.MyPrintStream;

import java.io.File;
import java.util.*;

public class Simulate {
	private double maxTrans, maxTorque;
	private Molecule[][][] lattice3;
//	PotentialLookup lookup;
	private double aConstant;
	private double totalEnergy;
	final static double scaleFactor = .25; 
	private double cutOff;
	private double torque;
	private double trans;
	private Random r;
	public final static boolean DEBUG = false;
	private Vector<PotentialLookup> potentials;
//	private MyPrintStream mpsLog;
//	private MyPrintStream mpsTransError;
	private boolean printPairwiseForcesToLog = false;
	/**
	 * 
	 * @param ljl
	 * @param cutOff distance cannot be less than this
	 * @param aConstant
	 */
	public Simulate(PotentialLookup lookup, double cutOff, double aConstant)
	{		
		this.aConstant = aConstant;
		totalEnergy = 0;
		
//		this.lookup = lookup;
		this.cutOff = cutOff;
		potentials = new Vector<PotentialLookup>();
		potentials.add(lookup);
//		mpsLog = new MyPrintStream(new File("simulate.log"));
//		File transLog = new File("transError.log");
//		int idx = 0;
//		if(transLog.exists())
//			do {
//				transLog = new File("transError--" + ++idx + ".log");
//			} while(transLog.exists());
//		mpsTransError = new MyPrintStream(transLog);
//		printInitialPotentials();
	}
//	private void printInitialPotentials() {
//		for(PotentialLookup pot : potentials) {
//			double x = 2;
//			double step = 0.01;
//			while(x < 10) {
//				x += step;
//				mpsLog.println(x + "\t" + pot.lookupForce(x) + "\t" + pot.lookupPotential(x));
//			}
//		}
//	}
	public void setMaxTrans(double maxTrans) { this.maxTrans = maxTrans; }
	public void setMaxTorque(double maxTorque) { this.maxTorque = maxTorque; }
	public void calculateTotalEnergy()
	{
		double energy = 0;
		
		for(int c = 0; c < getLattice3().length; c++)
			
			for(int b = 0; b < getLattice3()[c].length; b++)
				
				for(int a = 0; a < getLattice3()[c][b].length; a++)
					
					if( (a + b + c) % 2 == 0)
					{
						energy += calculateEnergyofFirstShell(getLattice3()[a][b][c]);
					}
		
		totalEnergy = energy / 2;
	}
	boolean testForContact(Molecule cur, Molecule surr, double cutoff) {
		for(int i = 0; i < cur.getLigands().length; i++)
			 for(int j = 0; j < surr.getLigands().length; j++)
				 if(JVector.subtract(cur.getLigands()[i].getPosition(), surr.getLigands()[j].getPosition()).length() < cutoff) { return true; }
				 
		return false;
	}
	public void resolve(Molecule cur, Molecule surr) {
		switch(new Random().nextInt(2))
		//switch(0)
		{
		case 0: 
			try {
				translate(surr, cur);
			} catch (Exception e) {
				copyOneToTwo(cur, surr);
				try {
					translate(surr, cur);
				} catch (Exception e1) {
				}
			}
			break;
		case 1: 
			try {
				torque(surr, cur);
			} catch (Exception e) {
				copyOneToTwo(cur, surr);
				try {
					torque(surr, cur);
				} catch (Exception e1) {
				}
			}
			break;
		}
	}
	public int removeContacts(double cutoff, int walkLength) {
		// initial variables
		
		Molecule cur = null;
		Molecule[] surr;
		int moves = 0, fails = 0, exitFromFails = 100, totalMoves = 0;
		int a, b, c, i;
		boolean hasContact;
		// after n consecutive failed tests raster through lattice testing for contacts
		while(fails < exitFromFails) {
			// select lattice point
			cur = null;
			do {
				a=r.nextInt(getLattice3().length);
				b=r.nextInt(getLattice3().length);
				c=r.nextInt(getLattice3().length);
				cur = getLattice3()[a][b][c];
			} while(cur == null);
			moves = 0;
			while(moves < walkLength) {
				// get surrounding
				surr = getSurrounding(cur).getSurroundingCBr4();
				
				hasContact = false;
				// test for br-br contacts < 3angstroms
				for(i = 0; i < surr.length; i++) {
					if(testForContact(cur, surr[i], cutoff)) {
						hasContact = true;
						break;
					}
				}
				// increment the fail counter if a contact does not exist
				if(!hasContact) {
					fails++;
					break;
				}
				// however if a contact does exist, reset the fail counter
				else { fails = 0; }
				// resolve conflict
				resolve(cur, surr[i]);
				// return molecules
				Molecule[] temp = {cur, surr[i]};
				returnMoleculesToOriginalPositions(temp);
				moves++;
				totalMoves++;
				// move to next lattice point while #moves < walkLength
				cur = surr[i];
			}
		}
		do {
			hasContact = false;
			for(a = 0; a < getLattice3().length; a++) {
				for(b = 0; b < getLattice3()[a].length; b++) {
					for(c = 0; c < getLattice3()[a][b].length; c++) {
						if((a+b+c) %2 == 0) {
							cur = getLattice3()[a][b][c];
							moves = 0;
							while(moves < walkLength) {
								// get surrounding
								surr = getSurrounding(cur).getSurroundingCBr4();
								
								hasContact = false;
								// test for br-br contacts < 3angstroms
								for(i = 0; i < surr.length; i++) {
									if(testForContact(cur, surr[i], cutoff)) {
										hasContact = true;
										break;
									}
								}
								// break the while loop if a contact does not exist
								if(!hasContact) { break; }
								else { hasContact = true; }
								// resolve conflict
								resolve(cur, surr[i]);
								moves++;
								totalMoves++;
								// move to next lattice point while #moves < walkLength
								cur = surr[i];
							}
						}
					}
				}
			}
		}
		while(hasContact);
		System.out.println("Total moves: " + totalMoves);
		return totalMoves;
	}
	public void MD(int numCycles) {
		for(int i = 0; i < numCycles; i++) {
			for(int a = 0; a < getLattice3().length; a++) {
				for(int b = 0; b < getLattice3()[a].length; b++) {
					for(int c = 0; c < getLattice3()[a][b].length; c++) {
						if(getLattice3()[a][b][c] != null) {
							relax(getLattice3()[a][b][c]);
							if(getLattice3()[a][b][c].getCenter().getPosition().i == Double.NaN || getLattice3()[a][b][c].getCenter().getPosition().j == Double.NaN || getLattice3()[a][b][c].getCenter().getPosition().k == Double.NaN) {
								System.out.println("dammit");
							}
						}
					}
				}
			}
		}

	}
	public void disorder() {
		JVector move = new JVector();
		JVector rot = new JVector();
		
		for(int a = 0; a < getLattice3().length; a++) {
			for(int b = 0; b < getLattice3()[a].length; b++) {
				for(int c = 0; c < getLattice3()[a][b].length; c++) {
					if((a+b+c)%2==0) {
						move.setI(r.nextDouble());
						move.setJ(r.nextDouble());
						move.setK(r.nextDouble());
						rot.setI(r.nextDouble());
						rot.setJ(r.nextDouble());
						rot.setK(r.nextDouble());
						try {
							move = JVector.multiply(move.unit(), maxTrans);
							rot = JVector.multiply(rot.unit(), maxTorque);
						} catch (Exception e) { e.printStackTrace(); }
						getLattice3()[a][b][c].translate(move);
						getLattice3()[a][b][c].rotate(rot, getLattice3()[a][b][c].getCenter().getPosition(), rot.length());
					}
				}	
			}	
		}
	}
	/**
	 * 
	 * @param numSteps
	 * @param numWalksPerMolecule
	 * @param enerCutoff interacting molecules must have an interaction potential greater than this value
	 * @param cutoff utilize the energy cutoff parameter
	 * @param disorderBetweenCycles randomly move the molecules between walk cycles
	 */
	public void walk(int numSteps, int numWalksPerMolecule, double enerCutoff, boolean cutoff, boolean disorderBetweenCycles) {
		int	max = getLattice3().length,
			a = 0, 
			b = 0, 
			c = 0, 
			idx = 0,
			totalCount = 0,
			count = 0;
		
		Molecule worst = null;

		double numMols = (int) (Math.pow(getLattice3().length/2, 3) * 4);
//		System.out.println("Target lattice energy: " + -1 * numMols * 0.763 + " eV");
		calculateTotalEnergy();
		System.out.println("Starting lattice energy: " + getTotalEnergy() + " eV");
		Vector<JVector> listOfPositions = new Vector<JVector>();
		for(a = 0; a < getLattice3().length; a++) {
			for(b = 0; b < getLattice3()[a].length; b++) {
				for(c = 0; c < getLattice3()[a][b].length; c++) {
					if(getLattice3()[a][b][c] != null) {
						listOfPositions.add(new JVector(a, b, c));
					}
				}
			}
		}
		for(int walkIdx = 0; walkIdx < numWalksPerMolecule; walkIdx++) {
			Vector<JVector> positions = new Vector<JVector>();
			for(int i = 0; i < listOfPositions.size(); i++) {
				positions.add((JVector) listOfPositions.get(i).clone());
			}
			idx = 0;
			if(disorderBetweenCycles) { disorder(); }				

			JVector pos;
			while(positions.size() > 0) {
				pos = positions.remove(r.nextInt(positions.size()));
				a = (int) pos.i;
				b = (int) pos.j;
				c = (int) pos.k;

				isNaN(getLattice3()[a][b][c].getCenter().getPosition().i);
				if(cutoff) { 
					worst = randomWalkHelper(getLattice3()[a][b][c], enerCutoff); 
				}
				else { 
					worst = randomWalkHelper(getLattice3()[a][b][c]); 
				}
				count = 1;
				while(count < numSteps) {
					count++;
					totalCount++;
					if(worst == null) { count = numSteps; }
					else {
						isNaN(worst.getCenter().getPosition().i);
						if(cutoff) { 
							worst = randomWalkHelper(worst, enerCutoff); 
						}
						else { 
							worst = randomWalkHelper(worst); 
						}
						isNaN(worst.getCenter().getPosition().i);
					} 
				}
			}
			calculateTotalEnergy();
			System.out.println("Total lattice energy after step: " + (walkIdx+1) + " of " + numWalksPerMolecule + ": " + getTotalEnergy() + " eV");
		}
		System.out.println("total walk steps: " + totalCount);

	}
	/**
	 * 
	 * @param walkLength
	 * @param numWalks
	 * @param enerCutoff
	 */
	public void walk(int walkLength, int numWalks, double enerCutoff)
	{
		int max = getLattice3().length;
		int idx = 0;
		for(int walkIdx = 0; walkIdx < numWalks; walkIdx++) {
			for(int i = 0; i < max; i++)
			{
				int a = 0;
				if(i == max-1) { a = max-1; }
				else { a = Integer.reverseBytes(i)%(max-1); }
				for(int j = 0; j < max; j++)
				{
					int b = 0;
					if(j == max-1) {b = max-1; }
					else { b = Integer.reverseBytes(j)%(max-1); }
					for(int k = 0; k < max; k++)
					{
						int c = 0;
						if(k == max-1) { c = max-1; }
						else { c = Integer.reverseBytes(k)%(max-1); }
						
						if((a + b + c) % 2 == 0)
						{
							Molecule worst = randomWalkHelper(getLattice3()[a][b][c], enerCutoff);
							
							int count = 1;
							
							while(count < walkLength)
							{
								//System.out.println(count);
								count++;
								if(worst != null) { 
									worst = randomWalkHelper(worst, enerCutoff); 
								} else { 
									count = walkLength; 
								}
							}
							idx++;
						}
					}
				}
			}
		}
	}
	private Molecule randomWalkHelper(Molecule center, double enerCutoff)
	{
		if(center == null) return null;
		int transOrRot = (new Random()).nextInt(2);
		
			Molecule worst = getLeastFavorableMolecule(center, enerCutoff);
			if(worst == null) return null;
			if(worst.equals(center))
				return null;
			else
			{
				switch(transOrRot)
				//switch(0)
				{
				case 0: 
					try {
						translate(worst, center);
					} catch (Exception e) {
						copyOneToTwo(center, worst);
						try {
							translate(worst, center);
						} catch (Exception e1) {
						}
					}
					break;
				case 1: 
					try {
						torque(worst, center);
					} catch (Exception e) {
						copyOneToTwo(center, worst);
						try {
							torque(worst, center);
						} catch (Exception e1) {
						}
					}
					break;
				}
	
				Molecule[] temp = {worst, center};
				returnMoleculesToOriginalPositions(temp);
				return worst;
			}
	}
//	private void writeToLog(Object obj) {
//		mpsLog.println(obj.toString());
//	}
	private Molecule randomWalkHelper(Molecule center) {
		if(center == null) return null;
		int transOrRot = r.nextInt(2);
		
			Molecule worst = getLeastFavorableMolecule(center);
			if(worst == null) return null;
			if(worst.equals(center))
				return null;
			else {
				switch(transOrRot) {
				case 0: 
					try {
						translate(worst, center);
						isNaN(worst.getCenter().getPosition().i);
						isNaN(center.getCenter().getPosition().i);
					} catch (Exception e) {
						copyOneToTwo(center, worst);
						try {
							System.out.println("trans err");
							translate(worst, center);
						} catch (Exception e1) {
						}
					}
					break;
				case 1: 
					try {
						torque(worst, center);
						isNaN(worst.getCenter().getPosition().i);
						isNaN(center.getCenter().getPosition().i);
					} catch (Exception e) {
						copyOneToTwo(center, worst);
						try {
							System.out.println("torque err");
							torque(worst, center);
						} catch (Exception e1) {
						}
					}
					break;
				}
	
				Molecule[] temp = {worst, center};
				returnMoleculesToOriginalPositions(temp);
				return worst;
			}
	}
	public void firstShellWalk(int walkLength, int numWalks) {
		int max = getLattice3().length;
		int idx = 0;
		for(int walkIdx = 0; walkIdx < numWalks; walkIdx++) {
			for(int i = 0; i < max; i++)
			{
				int a = 0;
				if(i == max-1) { a = max-1; }
				else { a = Integer.reverseBytes(i)%(max-1); }
				for(int j = 0; j < max; j++)
				{
					int b = 0;
					if(j == max-1) {b = max-1; }
					else { b = Integer.reverseBytes(j)%(max-1); }
					for(int k = 0; k < max; k++)
					{
						int c = 0;
						if(k == max-1) { c = max-1; }
						else { c = Integer.reverseBytes(k)%(max-1); }
						
						if((a + b + c) % 2 == 0)
						{
							Molecule worst = firstShellWalkHelper(getLattice3()[a][b][c], walkIdx);
							
							int count = 1;
							
							while(count < walkLength)
							{
								//System.out.println(count);
								count++;
								if(worst != null) { 
									worst = firstShellWalkHelper(worst, walkIdx); 
								} else { 
									count = walkLength; 
								}
							}
							idx++;
						}
					}
				}
			}
		}
	}
	private Molecule firstShellWalkHelper(Molecule center, int walkIdx)
	{
		if(center == null) return null;
		int transOrRot = (new Random()).nextInt(2);
		
			Molecule worst = getLeastFavorableMolecule(center, 0);
			if(worst == null) return null;
			if(worst.equals(center))
				return null;
			else
			{
				switch(transOrRot)
				{
				case 0: 
					try {
						translate2(center);
					} catch (Exception e) {
						return null;
//						copyOneToTwo(center, worst);
//						try {
//							translate2(center);
//						} catch (Exception e1) {
//						}
					}
					break;
				case 1: 
					try {
						torque2(center);
					} catch (Exception e) {
						copyOneToTwo(center, worst);
						try {
							torque2(center);
						} catch (Exception e1) {
						}
					}
					break;
				}
	
				Molecule[] temp = {worst, center};
				returnMoleculesToOriginalPositions(temp);
				return worst;
			}
	}
	private void copyOneToTwo(Molecule one, Molecule two)
	{
		JAtom[] ligands = one.getLigands();
		
		// one - two goes from two to one
		JVector difference = JVector.subtract(one.getCenter().getPosition(), two.getCenter().getPosition());
		
		two.translate(difference);
		
		two.setLigands(ligands);
		
		two.translate(JVector.multiply(difference, -1));
	}

	private void translate2(Molecule one) throws Exception
	{
		JVector[] forces = new JVector[13];
		JVector cur;
		Molecule[] surrounding = getSurrounding(one).getSurroundingCBr4();
		Molecule[] all = new Molecule[surrounding.length+1];
		all[0] = one;
		for(int i = 1; i < all.length; i++) {
			all[i] = surrounding[i-1];
		}
		// init all the force vectors
		for(int i = 0; i < all.length; i++) {
			forces[i] = new JVector();
		}
		for(int i = 0; i < all.length; i++) {
			for(int j = i; j < all.length; j++) {
				if(i==j) { continue; }
				cur = calculatePairwiseForceVector(all[i], all[j]);
				forces[i] = JVector.add(forces[i], JVector.multiply(cur, .5));
				forces[j] = JVector.add(forces[j], JVector.multiply(cur, -.5));
			}
		}
		
		// find the largest force
		double biggest = 0;
		double len = 0;
		for(int i = 0; i < all.length; i++) {
			len = forces[i].length();
			if(biggest < len) { biggest = len; }
		}
		for(int i = 0; i < all.length; i++) {
			forces[i] = JVector.multiply(forces[i], scaleFactor/biggest);
		}
		// translate all the molecules
		for(int i = 0; i < all.length; i++) {
			all[i].translate(forces[i]);
		}
	}	
	private void torque2(Molecule one) throws Exception
	{
		JVector[] forces = new JVector[13];
		JVector cur;
		Molecule[] surrounding = getSurrounding(one).getSurroundingCBr4();
		Molecule[] all = new Molecule[surrounding.length+1];
		all[0] = one;
		for(int i = 1; i < all.length; i++) {
			all[i] = surrounding[i-1];
		}
		// init all the force vectors
		for(int i = 0; i < all.length; i++) {
			forces[i] = new JVector();
		}
		for(int i = 0; i < all.length; i++) {
			for(int j = i; j < all.length; j++) {
				if(i==j) { continue; }
				cur = calculatePairwiseTorqueVector(all[i], all[j]);
				forces[i] = JVector.add(forces[i], JVector.multiply(cur, .5));
				forces[j] = JVector.add(forces[j], JVector.multiply(cur, -.5));
			}
		}
		
		// find the largest force
		double biggest = 0;
		double len = 0;
		for(int i = 0; i < all.length; i++) {
			len = forces[i].length();
			if(biggest < len) { biggest = len; }
		}
		for(int i = 0; i < all.length; i++) {
			forces[i] = JVector.multiply(forces[i], scaleFactor/biggest);
		}
		// rotate all the molecules
		for(int i = 0; i < all.length; i++) {
			all[i].rotate(forces[i], all[i].getCenter().getPosition(), forces[i].length());
		}
	}
	private void isNaN(double d) {
		if(Double.isNaN(d)) {
			System.out.println("Val is NaN.");
		}
	}
	private void translate(Molecule one, Molecule two) throws Exception
	{
		JVector force;
		try { 
			force = calculatePairwiseForceVector(one, two);
		}
		catch (Exception e1) { throw new Exception(); }
		isNaN(force.i);
		if(force.length() > maxTrans)
		{
			try { 
				force = JVector.multiply(force.unit(), maxTrans); 
			}
			catch (Exception e) {
				force = JVector.zero;
			}
		}
		force.multiply(.5);
		two.translate(force);

		force.multiply(-1);
		one.translate(force);

//		if(DEBUG) {
//			JVector interMolec = JVector.subtract(one.getCenter().position, two.getCenter().position);
//			double distBefore = interMolec.length();
//			System.out.println("Mol 1 center: " + one.getCenter());
//			System.out.println("Mol 2 center: " + two.getCenter());
//			writeToLog("Mol 1 center: " + one.getCenter());
//			writeToLog("Mol 2 center: " + two.getCenter());
//			interMolec = JVector.subtract(one.getCenter().position, two.getCenter().position);
//			double distAfter = interMolec.length();
//	
//			if(distBefore > distAfter) {
//				printPairwiseForcesToLog = true;
//				double potential = calculatePairwisePotential(one, two);
//				System.out.println("Distance shrunk.  Before: " + distBefore + "  After:" + distAfter + "\tPotential: " + potential);
//				mpsTransError.println("10\n");
//				mpsTransError.print(one.toStringForXYZ());
//				mpsTransError.print(two.toStringForXYZ());
//			} else
//				printPairwiseForcesToLog = false;
//		}
	}
	
	private JVector calculatePairwiseForceVector(Molecule one, Molecule two) throws Exception
	{
		JVector tmpForce, interAtom = null, interMolec = null, force = new JVector(0, 0, 0);
		JAtom[] atomsOnOne = one.getAtoms();
		JAtom[] atomsOnTwo = two.getAtoms();
		
//		if(printPairwiseForcesToLog) {
//			interMolec = JVector.subtract(one.getCenter().position, two.getCenter().position);
//		}
//		if(DEBUG) {
//			JVector ccDist = JVector.subtract(atomsOnOne[0].position, atomsOnTwo[0].position);
//			double ccDistVal = ccDist.length();
//			if(ccDistVal < 4) {
//				writeToLog("ccDist: " + ccDist.length() + "\t" + atomsOnOne[0].position.toString() + "\t" + atomsOnTwo[0].position.toString());
//			}
//			if(interMolec.length() < .1) {
//				return force;
//			}
//		}
		for(int j = 0; j < atomsOnOne.length; j++)
		{
			JAtom atomOne = atomsOnOne[j];
	
			for(int j1 = 1; j1 < 5; j1++)
			{
				JAtom atomTwo = atomsOnTwo[j1];
	
				try { 
					tmpForce = calculateForceVector(atomOne, atomTwo);
					force.add_inPlace(tmpForce);

//					if(printPairwiseForcesToLog) {
//						interAtom = JVector.subtract(atomOne.position, atomTwo.position);
//						mpsLog.println(atomOne + "\t" + atomTwo + "\t" + tmpForce + "\t" + interAtom + "\t" + 
//								interAtom.length() + "\t" + interMolec + "\t" + interMolec.length());
//					}
				}
				catch (Exception e) { throw new Exception(); }
				isNaN(force.i);
			}
		}
//		writeToLog(force);
		return force;
	}

	private JVector calculateForceVector(JAtom center, JAtom other) throws Exception
	{
		JVector direction = JVector.subtract(center.getPosition(), other.getPosition());
		
		double force;
		double length = direction.length();
		
		force = calculateForce(center.getZ(), other.getZ(), length); 

		
//		if(DEBUG) {
//			double potential = calculatePotential(center, other, length);
//			if(potential > 0 && force < 0) 
//				mpsLog.println("Potential: " + potential + "\tForce: " + force + "\tcenter: " + center + "\tother: " + other);
//		}
			
		direction = JVector.multiply(direction.unit(), force);
		direction.multiply(-1);
		return direction;
	}

	private void torque(Molecule one, Molecule two) throws Exception {
		JVector torque = null;
		try { torque = calculatePairwiseTorqueVector(one, two); }
		catch (Exception e1) { 
			try { throw new Exception(); }
			catch (Exception e) { throw new Exception(); }
		}
		double anglePos = torque.length();
		if(anglePos > maxTorque) {
			try { anglePos = JVector.multiply(torque.unit(), maxTorque).length(); }
			catch (Exception e) { e.printStackTrace(); }
		}
		
		if(anglePos != 0) {
			double angleNeg = -anglePos;
			one.rotate(torque, one.getCenter().getPosition(), anglePos);
			two.rotate(torque, two.getCenter().getPosition(), angleNeg);
		}
	}
	
	private JVector calculatePairwiseTorqueVector(Molecule one, Molecule two) throws Exception {
		JVector torque = new JVector(0, 0, 0);
		
		JAtom[] atomsOnOne = one.getAtoms();
		
		JAtom[] atomsOnTwo = two.getAtoms();
		
		for(int j = 0; j < atomsOnOne.length; j++)
	
			for(int j1 = 0; j1 < atomsOnTwo.length; j1++) {
	
				try { 
					torque = JVector.add(torque, calculateTorqueVector(one.getCenter(), atomsOnOne[j], atomsOnTwo[j1])); 
				}
				catch (Exception e) { 
					throw new Exception(); 
				}
			}
		
		return torque;
	}

	private JVector calculateTorqueVector(JAtom center, JAtom atomOne, JAtom atomTwo) throws Exception {
		JVector F = JVector.subtract(atomOne.getPosition(), atomTwo.getPosition());
		
		JVector r = JVector.subtract(atomOne.getPosition(), center.getPosition());
			
		try { F = F.unit(); }
		catch (Exception e) { e.printStackTrace(); }
		
		double force = 0;
		try { 
			force = calculateForce(atomOne.getZ(), atomTwo.getZ(), F.length());
			if(force == 0)
				return new JVector(0, 0, 0);
		}
		catch (Exception e) { throw new Exception(); }
		
		F = JVector.multiply(F.unit(), force);
		
		JVector torque = JVector.cross(r, F);
		
		return torque;
	}

	private double calculateForce(int Z1, int Z2, double distance) throws Exception {
//		double distance = JVector.subtract(center.position, other.position).length();
		
		if(distance < cutOff)
			throw new Exception();
		for(PotentialLookup potLook : potentials)
			if(potLook.check(Z1, Z2)) {
				double forceVal = potLook.lookupForce(distance);
				return forceVal;
			}
		
		return 0;
	}

	private double calculatePotential(JAtom center, JAtom other, double distance) throws Exception {
//		double distance = JVector.subtract(center.position, other.position).length();
		
		if(distance < cutOff)
			throw new Exception();
		for(PotentialLookup potLook : potentials)
			if(potLook.check(center.getZ(), other.getZ())) {
				double potVal = potLook.lookupPotential(distance);
				return potVal;
			}
		
		return 0;
	}

	/**
	 * Method to determine the face that the molecule is closest to and if the molecule is
	 * "on the wrong side of the world" then subtract the size of the box from the coordinate
	 * By "on the wrong side of the world" I mean that if the molecule has int coordinates of 
	 * (0, 0, 0) but has been translated by (-.01, -.01, -.01) then the molecule will really be at
	 * (if the box is 44.1 angstroms on a side) (44.09, 44.09, 44.09).  For the purposes of this, I
	 * will need to subtract the box size from the molecule to put it back at (-.01, -.01, -.01).
	 * 
	 * @param central The central tetrahedron
	 * @param axis the coordinate switch: 0 - x, 1 - y, 2 - z
	 * @return	The face that the molecule is closest to.  2 - not close to a face, 0 - min face, 1 - max face
	 */
	private int determineFace(Molecule central, int axis) {
		// get the parameters from the lattice3
		int aUnits = getLattice3().length;
		
		int face = 2;	// define the face parameter, 2 = not near a face
		
		double coordinate_d = 0;	// define the real coordinate of the molecule
		
		int coordinate_i = 0;	// define the lattice3 position of the molecule

		// switch block to pick the correct coordinate
		switch(axis) {
		case 0: coordinate_d = central.getCenter().getPosition().getI(); break;
		case 1: coordinate_d = central.getCenter().getPosition().getJ(); break;
		case 2: coordinate_d = central.getCenter().getPosition().getK(); break;
		}
		
		// take the double coordinate of the molecule and convert it into lattice3 coordinates
		// i.e. coordinate_d = 8.82, coordinate_i = 2
		coordinate_i = (int) Math.round(coordinate_d * 2 / aConstant) % aUnits;

		// test to see if the integer coordinate is near the min face
		if(coordinate_i <= 1) {
			face = 0;
			// test to see if the coordinate is on the "wrong side of the world" and if so set the 
			// boolean subtract flag to true
			if(testForWorldChanger(coordinate_i, coordinate_d))
				changeSidesOfWorld(central, axis);
		}
		// test to see if the integer coordinate is near the max face
		else if(coordinate_i >= aUnits - 1)
			face = 1;
		
		// return the face which the molecule is close to
		return face;
	}
	
	private boolean testForWorldChanger(int latticePos, double realPos) {
		if(latticePos == 0 && realPos > aConstant)
			return true;

		return false;
	}
	
	
	/**
	 * Method to move a molecule to the "correct side of the world"
	 * By "on the wrong side of the world" I mean that if the molecule has int coordinates of 
	 * (0, 0, 0) but has been translated by (-.01, -.01, -.01) then the molecule will really be at
	 * (if the box is 44.1 angstroms on a side) (44.09, 44.09, 44.09).  For the purposes of this, I
	 * will need to subtract the box size from the molecule to put it back at (-.01, -.01, -.01).
	 * @param worldChanger	The molecule that is to be moved
	 * @param axis	The axis that the molecule is close to.  0 - x, 1 - y, 2 - z
	 */
	private void changeSidesOfWorld(Molecule worldChanger, int axis)
	{
		// get the parameters from the lattice3
		int aUnits = getLattice3().length;
		
		JVector translate = new JVector(0, 0, 0);	// define a new vector object
		
		// switch block to decide which direction the molecule needs to translate
		switch(axis)
		{
		case 0: translate.setI(-aConstant * aUnits / 2); break; // x
		case 1: translate.setJ(-aConstant * aUnits / 2); break;	// y
		case 2: translate.setK(-aConstant * aUnits / 2); break;	// z
		}
		// translate the molecule in the correct direction
		worldChanger.translate(translate);
	}
	
	private double determineTranslateAmount(int face, int position)
	{
		int aUnits = getLattice3().length;
		
		double toTranslate = 0;
		
		switch(face)
		{
		case 0:	
			if(position == aUnits - 1)
				toTranslate = - aConstant * aUnits / 2;
			break;
		case 1:
			if(position == 0)
				toTranslate = + aConstant * aUnits / 2;
			break;
		}
		
		return toTranslate;
	}
	
	private double calculatePairwisePotential(Molecule one, Molecule two)
	{
		double potential = 0;
		
		JAtom[] atomsOnOne = one.getAtoms();
		
		JAtom[] atomsOnTwo = two.getAtoms();
		
		for(int j = 0; j < 5; j++)
		{
			JAtom atomOne = atomsOnOne[j];

			for(int j1 = 0; j1 < 5; j1++)
			{
				JAtom atomTwo = atomsOnTwo[j1];
				
				potential += calculateEnergy(atomOne, atomTwo);
			}
		}
		
		return potential;
	}
	
	/**
	 * Method that looks at the first shell of the molecule specified by the 'center' variable and calculates
	 * the energy of interaction between the central molecule and the surrounding molecules.  The molecule with
	 * the most unfavorable energy of interaction is then returned.
	 * @param center	The central molecule
	 * @return	The molecule with the most unfavorable energy of interaction
	 */
	Molecule getLeastFavorableMolecule(Molecule center, double enerCutoff)
	{
		DoubleLinkedListMolecule listCBr4 = getSurrounding(center);
		//if(listCBr4.length < 12) return null;
		//System.out.println("length of double linked list: " + listCBr4.length);
		//for(int i = 0; i < listCBr4.length; i++) {
		//	System.out.println(listCBr4.removeHead().getValue());
		//}
		//System.out.println("length of double linked list: " + listCBr4.length);
		// remove the molecule with the worst energy of interaction
		//TwoNodeCBr4 worstNode = listCBr4.removeLast();
		
		
		// get the linked list as an array
		Molecule[] temp = new Molecule[1];
		// return the molecules to their original positions
		Molecule worst = null;
		double ener = 0;
		Molecule[] list = null;
		int Z1, Z2;
		JAtom[] centerAtoms, worstAtoms;
		while(listCBr4.length != 0) {
			worst = listCBr4.removeHead().value;
			centerAtoms = center.getAtoms();
			worstAtoms = worst.getAtoms();
			for(int i = 0; i < centerAtoms.length; i++) {
				Z1 = centerAtoms[i].getZ();
				for(int j = 0; j < worstAtoms.length; j++) {
					Z2 = centerAtoms[j].getZ();
					for(PotentialLookup potLook : potentials) {
						if(!potLook.check(Z1, Z2))
							continue;
						
						ener = potLook.lookupPotential(JVector.subtract(center.getLigands()[i].getPosition(), worst.getLigands()[j].getPosition()).length());
						if(ener >= enerCutoff) { 
							list = listCBr4.getSurroundingCBr4();
							returnMoleculesToOriginalPositions(list);
							return worst; 
						}
					}
				}
			}
			temp[0] = worst;
			returnMoleculesToOriginalPositions(temp);
		}
		
		return null;
	}
	Molecule getLeastFavorableMolecule(Molecule center )
	{
		DoubleLinkedListMolecule listCBr4 = getSurrounding(center);
		//if(listCBr4.length < 12) return null;
		//System.out.println("length of double linked list: " + listCBr4.length);
		//for(int i = 0; i < listCBr4.length; i++) {
		//	System.out.println(listCBr4.removeHead().getValue());
		//}
		//System.out.println("length of double linked list: " + listCBr4.length);
		// remove the molecule with the worst energy of interaction
		TwoNodeMolecule worstNode = listCBr4.removeLast();
		while(JVector.distance(center.getCenter().getPosition(), worstNode.value.getCenter().getPosition()) < .1) {
			listCBr4.insertHead(worstNode);
			worstNode = listCBr4.removeLast();
		}
		
		// get the linked list as an array
		Molecule[] list = listCBr4.getSurroundingCBr4();
		//Molecule[] temp = new Molecule[1];
		// return the molecules to their original positions
		returnMoleculesToOriginalPositions(list);

		return worstNode.value;
	}
	public DoubleLinkedListMolecule getSurrounding(Molecule center)
	{
		//Sorted Linked List to hold the thirteen molecules
		DoubleLinkedListMolecule listCBr4 = new DoubleLinkedListMolecule(13);
		
		// The lattice position of the molecule
		JVector molecularCenter_i = JVector.multiply(center.getCenter().getPosition(), 2 / aConstant).roundInt();
		
		// declare the faces array that will hold the faces which the molecule is closest to
		int[] faces = new int[3];
		//  the array 
		for(int i = 0; i < faces.length; i++)
			faces[i] = determineFace(center, i);

		int aMax = (int) molecularCenter_i.getI() + 1;
		int bMax = (int) molecularCenter_i.getJ() + 1;
		int cMax = (int) molecularCenter_i.getK() + 1;

		int x = 0;
		int y = 0;
		int z = 0;
		
		for(int c = cMax - 2; c <= cMax; c++)
		{
			// check to make sure the c index is inside the box
			z = checkIntCoordinate(c);
			for(int b = bMax - 2; b <= bMax; b++)
			{
				// check to make sure the b index is inside the box
				y = checkIntCoordinate(b);
				for(int a = aMax - 2; a <= aMax; a++)
				{
					// test to see if there is a molecule at the specified location
					//if( (a + b + c) % 2 == 0)
					{
						// test to see if the molecule specified by (a, b, c) is the same as the central molecule
						if(a == aMax-1 && b == bMax-1 && c == cMax-1)
							continue;
						// check to make sure the a index is inside the box
						x = checkIntCoordinate(a);

						if(getLattice3()[x][y][z] == null) { continue; }
						// get the molecule at the lattice location
						Molecule surrounding = getLattice3()[x][y][z];
						if(surrounding == null) { continue; }
						if(surrounding == center) { continue; }
						//System.out.println(surrounding);
						// get the center of the molecule
						JVector surroundingCenter = surrounding.getCenter().getPosition();
						
						// move the molecule to the correct side of the world if necessary
						if(x == 0 && testForWorldChanger(x, surroundingCenter.getI()))
							changeSidesOfWorld(surrounding, 0);
						if(y == 0 && testForWorldChanger(y, surroundingCenter.getJ()))
							changeSidesOfWorld(surrounding, 1);
						if(z == 0 && testForWorldChanger(z, surroundingCenter.getK()))
							changeSidesOfWorld(surrounding, 2);
						
						JVector translate = new JVector(0, 0, 0);
					
						// determine the amount that the molecule needs to translate
						translate.setI(determineTranslateAmount(faces[0], x));	// x
						translate.setJ(determineTranslateAmount(faces[1], y));	// y
						translate.setK(determineTranslateAmount(faces[2], z));	// z
						
						// translate the molecule by the amount defined above
						surrounding.translate(translate);
						
						// get the new center of the molecule
						surroundingCenter = surrounding.getCenter().getPosition();
						
						// calculate the energy of interaction between the selected molecule and the central molecule
						double energy = calculatePairwisePotential(center, surrounding);
						
						if(Double.isInfinite(energy) || energy > 10000) {
//							System.out.println("Infinite Energy");
							surrounding = null;
							surrounding = Lattice.fillInSpaces(lattice3, aConstant);
							energy = calculatePairwisePotential(center, surrounding);
						}
						
						// insert a new node into the sorted linked list
						listCBr4.insert(new TwoNodeMolecule(surrounding, energy, null, null));
					}
				}
			}
		}
		// return the molecule that is interacting most poorly with the central molecule
		return listCBr4;
	}
	
	/**
	 * Method to return molecules that have been moved to their original positions
	 * @param list	The list of molecules that are to be moved
	 */
	private void returnMoleculesToOriginalPositions(Molecule[] list)
	{
		for(int i = 0; i < list.length; i++)
		{
			// get the current position of the molecule
			JVector oldPosition = list[i].getCenter().getPosition();
			
			JVector newPosition = new JVector(0, 0, 0);

			// check each coordinate of the current position of the molecule and set it to the proper coordinate if
			// it is outside the boundaries of the box
			newPosition.setI(checkDoubleCoordinate(oldPosition.getI()));	// x
			newPosition.setJ(checkDoubleCoordinate(oldPosition.getJ()));	// y
			newPosition.setK(checkDoubleCoordinate(oldPosition.getK()));	// z
			
			// calculate the vector from oldPosition -> newPosition
			newPosition = JVector.subtract(newPosition, oldPosition);
			
			// translate the molecule by that amount
			list[i].translate(newPosition);
		}
	}
	
	
	public double calculateEnergyofFirstShell(Molecule center)
	{
		double energy = 0;
		
		DoubleLinkedListMolecule list = getSurrounding(center);
		
		TwoNodeMolecule[] removed = new TwoNodeMolecule[list.length];
		
		Molecule[] removedMolecules = new Molecule[list.length];
		
		int counter = 0;
		
		int listLength = list.length;
		
		while(counter < listLength)
		{
			removed[counter] = list.removeHead();
			
			removedMolecules[counter] = removed[counter].value;
			
			energy += removed[counter].key;
			
			counter++;
		}
		
		returnMoleculesToOriginalPositions(removedMolecules);
		
		return energy;
	}
	
	private double calculateEnergy(JAtom center, JAtom other)
	{
		for(PotentialLookup potLook : potentials) {
			if(potLook.check(center.Z, other.Z)) {
				double distance = JVector.subtract(center.getPosition(), other.getPosition()).length();
				return potLook.lookupPotential(distance);
			}
		}
		return 0;
	}

	private int checkIntCoordinate(int coordinate)
	{
		int a = getLattice3().length;
		
		return ((coordinate % a) + a) % a;
	}
	
	private double checkDoubleCoordinate(double coordinate)
	{
		double a = getLattice3().length * aConstant / 2;
		
		return ((coordinate % a) + a) % a;
	}

	public double getTotalEnergy() { return totalEnergy; }
	
	/**
	 * 
	 * @param whichType 0-translate, 1-rotate, 2-translate-rotate, 3-rotate-translate, 4-randomly translate or rotate
	 * 5-randomly trans, rot, trans-rot, rot-trans
	 */
	public void iterate(int numCycles, int whichType)
	{
		int max = getLattice3().length;
		for(int cycles = 0; cycles < numCycles; cycles++) {
			for(int i = 0; i < max; i++)
			{
				int a = 0;
				if(i == max-1) { a = max-1; }
				else { a = Integer.reverseBytes(i)%(max-1); }
				for(int j = 0; j < max; j++)
				{
					int b = 0;
					if(j == max-1) {b = max-1; }
					else { b = Integer.reverseBytes(j)%(max-1); }
					for(int k = 0; k < max; k++)
					{
						int c = 0;
						if(k == max-1) { a = max-1; }
						else { a = Integer.reverseBytes(k)%(max-1); }
						
						if((a + b + c) % 2 == 0)
						{
							switch(whichType)
							{
							case 0:
								stage2VectorField(getLattice3()[a][b][c]);
								break;
							case 1:
								stage3Torque(getLattice3()[a][b][c]);
								break;
							case 2:
								stage2VectorField(getLattice3()[a][b][c]);
								stage3Torque(getLattice3()[a][b][c]);
								break;
							case 3:
								stage3Torque(getLattice3()[a][b][c]);
								stage2VectorField(getLattice3()[a][b][c]);
								break;
							case 4:
								switch(new Random().nextInt(2))
								{
								case 0:
									stage2VectorField(getLattice3()[a][b][c]);
									break;
								case 1:	
									stage3Torque(getLattice3()[a][b][c]);
									break;
								}
							case 5:
								switch(new Random().nextInt(4))
								{
								case 0:
									stage2VectorField(getLattice3()[a][b][c]);
									break;
								case 1:
									stage3Torque(getLattice3()[a][b][c]);
									break;
								case 2:
									stage2VectorField(getLattice3()[a][b][c]);
									stage3Torque(getLattice3()[a][b][c]);
									break;
								case 3:
									stage3Torque(getLattice3()[a][b][c]);
									stage2VectorField(getLattice3()[a][b][c]);
									break;
								}
							case 6: relax(getLattice3()[a][b][c]);
							}
							
						}
					}
				}
			}
		}
	}
	
	public void stage2VectorField(Molecule molecule)
	{
		DoubleLinkedListMolecule center = getSurrounding(molecule);
		Molecule[] thirteen = center.getSurroundingCBr4();
		JVector force = calculateForceVectorFor13(molecule, thirteen);
		
		if(force.length() > .5)
			try {
				force = JVector.multiply(force.unit(), 0.5);
			} catch (Exception e) {
				e.printStackTrace();
			}
		
		molecule.translate(force);
		returnMoleculesToOriginalPositions(thirteen);
	}
	
	public void stage3Torque(Molecule molecule)
	{
		DoubleLinkedListMolecule center = getSurrounding(molecule);
		Molecule[] thirteen = center.getSurroundingCBr4();
		JVector torque = calculateForceVectorFor13(molecule, thirteen);
		
		if(torque.length() > .5)
			try {
				torque = torque.unit();
			} catch (Exception e) {
				e.printStackTrace();
			}
		
		molecule.rotate(torque, molecule.getCenter().getPosition(), torque.length());

		returnMoleculesToOriginalPositions(thirteen);
	}
	public void relax(Molecule molecule) {
		if(molecule == null) {
			System.out.println("wtf");
		}
		Molecule[] surrounding = getSurrounding(molecule).getSurroundingCBr4();
		JVector[][] forces = calcTorqueAndTransFor13(molecule, surrounding);
		//double enerBefore = calculateEnergyofFirstShell(molecule);
		double maxLen = 0;
		for(int i = 0; i < forces.length; i++) {
			for(int j = 0; j < forces[i].length; j++) {
					if(forces[i][j].length() > maxLen) { maxLen = forces[i][j].length(); }
			}
		}
		double scalar = maxLen / .25;
		for(int i = 0; i < forces.length; i++) {
			for(int j = 0; j < forces[i].length; j++) {
				if(forces[i][j].i == Double.NaN ||forces[i][j].j == Double.NaN ||forces[i][j].k == Double.NaN) {
					System.out.println("dammit");
				}
					forces[i][j] = JVector.multiply(forces[i][j], 1./scalar);
			}
		}
		molecule.translate(forces[0][0]);
		molecule.rotate(forces[1][0], molecule.getCenter().getPosition(), forces[1][0].length());
		/*for(int j = 1; j < forces[0].length; j++) {
			surrounding[j-1].translate(forces[0][j]);
			surrounding[j-1].rotate(forces[1][j], surrounding[j-1].getCenter().position, forces[1][j].length());
		}
		double enerAfter = calculateEnergyofFirstShell(molecule);
		System.out.println("Energy lowered by: " + (enerBefore - enerAfter));*/
		returnMoleculesToOriginalPositions(surrounding);
	}
	/**
	 * 
	 * @param molecule
	 * @param surrounding
	 * @return 2d array:
	 * <br>[0] are translation vectors,
	 * <br>[1] are torque vectors
	 * <br>[][0] is the central molecule
	 * <br>[][1-12] are the surrounding molecules
	 */
	public JVector[][] calcTorqueAndTransFor13(Molecule molecule, Molecule[] surrounding) {
		JVector f, t;
		// [0] are translation vectors
		// [1] are torque vectors
		// [][0] is the central molecule
		JVector[][] tnt = new JVector[2][surrounding.length+1];
		tnt[0][0] = new JVector();
		tnt[1][0] = new JVector();
		for(int i = 1; i < surrounding.length+1; i++) {
			try {
				f = JVector.multiply(calculatePairwiseForceVector(molecule, surrounding[i-1]), .5);
				t = JVector.multiply(calculatePairwiseTorqueVector(molecule, surrounding[i-1]), .5);
				tnt[0][0] = JVector.add(tnt[0][0], f);
				tnt[1][0] = JVector.add(tnt[1][0], t);
				
				tnt[0][i] = JVector.multiply(f, -1);
				tnt[1][i] = JVector.multiply(t, -1);
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return tnt;
	}
	public JVector calculateForceVectorFor13(Molecule molecule, Molecule[] surrounding)
	{
		JVector force = new JVector();
		
		for(int i = 0; i < 12; i++)
		{
			try {
				force = JVector.add(force, calculatePairwiseForceVector(molecule, surrounding[i]));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		return force;
	}
	
	public JVector calculateTorqueVectorFor13(Molecule molecule, Molecule[] surrounding)
	{
		JVector torque = new JVector();
		
		for(int i = 0; i < 12; i++)
		{
			try {
				torque = JVector.add(torque, calculatePairwiseTorqueVector(molecule, surrounding[i]));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		return torque;
	}
	
	public int[] monteCarloOrientation(int numTestsPerMolecule) {
		double enerBefore, enerAfter;
		
		int numTotalUnphysicalContacts = 0;
		int numUnphysicalRemoved = 0;
		int a, b, c;
		Molecule molBefore;
		for(int i = 0; i < numTestsPerMolecule; i++) {
			Vector<JVector> positions = new Vector<JVector>();
			for(a = 0; a < lattice3.length; a++) {
				for(b = 0; b < lattice3[a].length; b++) {
					for(c = 0; c < lattice3[a][b].length; c++) {
						if(lattice3[a][b][c] != null) {
							positions.add(new JVector(a, b, c));
						}
					}
				}
			}
			
			numTotalUnphysicalContacts = 0;
			numUnphysicalRemoved = 0;
			
			
			while(positions.size() > 0) {
				JVector curPos = positions.remove(r.nextInt(positions.size()));
				a = (int) curPos.i;
				b = (int) curPos.j;
				c = (int) curPos.k;
				numTotalUnphysicalContacts++;
				molBefore = lattice3[a][b][c];
				do {
					enerBefore = calculateEnergyofFirstShell(molBefore);
					newOrientation(a, b, c);
					enerAfter = calculateEnergyofFirstShell(lattice3[a][b][c]);
				} while(enerAfter == enerBefore);
				if(enerAfter > enerBefore) {
					lattice3[a][b][c] = molBefore;
				} else {
					numUnphysicalRemoved++;
				}
			}
		}
		return new int[] {numTotalUnphysicalContacts, numUnphysicalRemoved};
	}
	private void newOrientation(int a, int b, int c) {
		JVector move = lattice3[a][b][c].getCenter().getPosition();
		
		int newChoice = r.nextInt(6);
		
		lattice3[a][b][c] = Lattice.getOrientation(newChoice);
		
		lattice3[a][b][c].translate(move);
	}
	public void monteCarloShift(int numShift) {
		
	}
	
	/**
	 * 
	 * @return 2d array of molecules.  First index corresponds to the 1d Molecule[] array from the to1dArray() call.
	 * <br>Second index are the twelve molecules around the central molecule.
	 */
	public Molecule[][] getSurrounding() {
		int numMolecules = 0;
		for(int i = 0; i < lattice3.length; i++)
			for(int j = 0; j < lattice3[i].length; j++)
				for(int k = 0; k < lattice3[i][j].length; k++)
					if(lattice3[i][j][k] != null) {
						numMolecules++;
					}
//		System.out.println("numMolecules: " + numMolecules);

		Molecule[][] surrounding = new Molecule[numMolecules][];
		int idx = 0;
		for(int i = 0; i < lattice3.length; i++) {
			for(int j = 0; j < lattice3[i].length; j++)
				for(int k = 0; k < lattice3[i][j].length; k++)
					if(lattice3[i][j][k] != null) {
						surrounding[idx++] = getSurrounding(lattice3[i][j][k]).getSurroundingCBr4();
//						System.out.println(idx + " of " + numMolecules);
					}
//			if(idx%1000 == 0)
//				System.out.println(idx + " of " + numMolecules);
		}
		
		return surrounding;
	}
	
	/**
	 * 
	 * @return 1d list of molecules in the lattice
	 */
	public Molecule[] to1dArray() {
		Vector<Molecule> mols = new Vector<Molecule>();
		for(Molecule[][] a : lattice3)
			for(Molecule[] b : a)
				for(Molecule c : b)
					if(c != null)
						mols.add(c);
		
		return mols.toArray(new Molecule[mols.size()]);
	}
	public void setRandom(Random r) { this.r = r; }
	public double getTorque() { return torque; }
	public void setTorque(double torque) { this.torque = torque; }
	public double getTrans() { return trans; }
	public void setTrans(double trans) { this.trans = trans; }
	public Molecule[][][] getLattice3() { return lattice3; }
	public void setLattice3(Molecule[][][] lattice3) { this.lattice3 = lattice3; }
	public Vector<PotentialLookup> getPotentialsLookup() {
		return potentials;
	}
	public void setPotentialsLookup(Vector<PotentialLookup> potentials) {
		this.potentials = potentials;
	}
	public void addPotentialLookup(PotentialLookup potLook) { 
		if(potentials == null || potentials.get(0) == null)
			potentials = new Vector<PotentialLookup>();
		
		potentials.add(potLook); 
	}
}
