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

import io.MyPrintStream;

import java.io.File;
import java.util.Vector;

import chemistry.JAtomTools;

public class PotentialTest {

	private Vector<Potential> potentials;
	private double timeStep = 0.0001;
	
	public JAtom[][] buildSurroundings(JAtom[] atoms) {
		Vector<JAtom[]> allSurroundings = new Vector<JAtom[]>();
		double maxBrBrDist = 7;
		double maxCBrDist = 2.5;
		
		for(JAtom atom1 : atoms) {
			Vector<JAtom> surroundings = new Vector<JAtom>();
			for(JAtom atom2 : atoms) {
				if(atom1 == atom2)
					continue;
				
				if((atom1.getZ() == 6 && atom2.getZ() == 35 ||
						atom1.getZ() == 35 && atom2.getZ() == 6) &&
						JVector.distance(atom1.getPosition(), atom2.getPosition()) < maxCBrDist)
					surroundings.add(atom2);
				
				if((atom1.getZ() == 35 && atom2.getZ() == 35) &&
						JVector.distance(atom1.getPosition(), atom2.getPosition()) < maxBrBrDist)
					surroundings.add(atom2);
			}
			allSurroundings.add(surroundings.toArray(new JAtom[surroundings.size()]));
		}
		
		return allSurroundings.toArray(new JAtom[allSurroundings.size()][]);
	}
	private void calcForces(JAtom[] atoms, JAtom[][] surroundings, JVector[] forces) {
		int Z1, Z2;
		JVector forceVec, pos1, pos2;
		double force;
		for(int i = 0; i < atoms.length; i++) {
			Z1 = atoms[i].getZ();
			pos1 = atoms[i].getPosition();
			for(int j = 0; j < surroundings[i].length; j++) {
				Z2 = surroundings[i][j].getZ();
				pos2 = surroundings[i][j].getPosition();
				for(Potential pot : potentials) {
					if(pot.check(Z1, Z2)) {
						forceVec = JVector.subtract(pos1, pos2);
						force = pot.calcF(forceVec.length());
						forceVec = JVector.multiply(forceVec.unit(), force);
						forces[i].add_inPlace(forceVec);
					}
				}
			}
		}
	}
	private void moveAtoms(JAtom[] atoms, JVector[] forces) {
		JVector force;
		double mass;
		double forceMagnitude;
		for(int i = 0; i < atoms.length; i++) {
			force = forces[i];
			mass = JAtomTools.getMass(atoms[i].getZ());
			forceMagnitude = 0.5 * mass * timeStep * timeStep;
			force = JVector.multiply(force, forceMagnitude);
			System.out.println("Amount to move: " + force.length());
			atoms[i].getPosition().add_inPlace(force);
		}
	}
	private void printPositions(MyPrintStream mps, JAtom[] atoms) {
		for(JAtom atom : atoms) {
			mps.println(atom.getZ() + "\t" + atom.getPosition().toTabString());
		}
	}
	public void run() {
		/* ************************** */
		/* HOOKE POTENTIAL PARAMETERS */
		/* ************************** */
		PotentialLookup hookeLook;
		double tolerance = 1e-4;
		
		HookePotential hooke;
		double cBrDist = 1.91;
		// nu = 267 cm-1 from handbook of chemistry and physics
		// mu = 12.01 * 79.9 / (12.01 + 79.9) / 1000 / 6.022e23 = 1.73e-26
		// k = (2*PI*3e10*267)^2 * 1.73e-26 = 43.82
		double cBrSpringConstant;
		cBrSpringConstant = 43.82;	// upper Limit
//		cBrSpringConstant /= 4;		// lower limit
		
		/* ********************************** */
		/* LENNARD JONES POTENTIAL PARAMETERS */
		/* ********************************** */
		PotentialLookup ljl;
		LennardJonesPotential ljp;
		double lookupPrecision = 0.001;
		int Z1 = 35;
		int Z2 = 35;
		double Z1Z2_distance = 3.8;
		double zeroPotentialDistance;
		double wellDepth = 0.034891;
		
		zeroPotentialDistance = Z1Z2_distance / Math.pow(2., (1./6.));
		ljp = new LennardJonesPotential(Z1, Z2, zeroPotentialDistance, wellDepth);
		hooke = new HookePotential(6, 35, cBrDist, cBrSpringConstant);

		JAtom center = new JAtom(6, new JVector(0, 0, 0));
		JVector ligandPos = new JVector(1, 1, 1).unit();
		ligandPos.multiply(1.91);
		JAtom ligand = new JAtom(35, ligandPos);
		IdealTetrahedron test = new IdealTetrahedron(center, ligand);
		
		Vector<JAtom> tempAtoms = new Vector<JAtom>();
		for(JAtom atom : test.getAtoms())
			tempAtoms.add(atom);
		
		JAtom[] atoms = new JAtom[tempAtoms.size()]; 
		atoms = tempAtoms.toArray(atoms);
//		while(true) {
		
		JAtom br = atoms[1];
		JVector brBrforce = new JVector();
		JVector betweenBr;
		double distBetweenBr;
		for(int i = 2; i < atoms.length; i++) {
			betweenBr = JVector.subtract(br.getPosition(), atoms[i].getPosition());
			distBetweenBr = betweenBr.length();
			brBrforce.add_inPlace(JVector.multiply(betweenBr.unit(), ljp.calcF(distBetweenBr)));
			
		}
		double brBrDistBefore = JVector.distance(atoms[1].getPosition(), atoms[2].getPosition());
		double brBrForceBefore = ljp.calcF(brBrDistBefore);
		double cBrDistBefore = JVector.distance(atoms[0].getPosition(), atoms[1].getPosition());
		double newRestingPosition = cBrDist-(brBrforce.length() / cBrSpringConstant);
		hooke = new HookePotential(6, 35, newRestingPosition, cBrSpringConstant);
		double cBrForceBefore = hooke.calcF(cBrDistBefore);
		
		
		JAtom[][] allSurroundings = buildSurroundings(atoms);
		for(int i = 0; i < allSurroundings.length; i++) {
			System.out.println("Atom type/idx: " + atoms[i].getZ() + "/" + atoms[i].atomID);
			for(int j = 0; j < allSurroundings[i].length; j++) {
				System.out.println("\tSurrounding atom type/idx: " + allSurroundings[i][j].getZ() + 
						"/" + allSurroundings[i][j].atomID);
			}
		}
		JVector[] force = new JVector[atoms.length];
		for(int i = 0; i < force.length; i++)
			force[i] = new JVector();
		potentials = new Vector<Potential>();
		potentials.add(ljp);
		potentials.add(hooke);
		MyPrintStream movie = new MyPrintStream(new File("test.xyz"));
		movie.println(atoms.length);
		movie.println();
		
		int idx = 0;
		while(idx < 100) {
			for(int i = 0; i < force.length; i++) {
				force[i].setI(0);
				force[i].setJ(0);
				force[i].setK(0);
			}
			calcForces(atoms, allSurroundings, force);
			moveAtoms(atoms, force);
			printPositions(movie, atoms);
//			movie.println();
			idx++;
		}
		movie.flush();
		movie.close();
//		}
//		System.out.println("Br-Br intramolecular force: " + ljp.calcF(3.118993));
//		System.out.println("C-Br intramolecular force: " + ljp.calcF(3.118993));

		double step = .05;
		double startVal = 1;
		double finishVal = 3;
		int numSteps = (int) Math.rint((finishVal - startVal) / step);
		JAtom c = atoms[0];
		br = atoms[1];
		JVector cBrVec = JVector.subtract(br.getPosition(), c.getPosition()),
				intraMolecVec;
		double potVal, forceVal;
		JAtom a1, a2;
		System.out.println("Distance\tPotentialEnergy\tForce:");
		for(int i = 0; i < numSteps; i++) {
			cBrVec = cBrVec.unit();
			cBrVec.multiply(startVal + i*step);
			br.setNewPos(cBrVec);
			forceVal = potVal = 0;
			for(int a = 0; a < atoms.length; a++) {
				a1 = atoms[a];
				if(a1 == br)
					continue;
				
				intraMolecVec = JVector.subtract(br.getPosition(), a1.getPosition());
				for(Potential pot : potentials)
					if(pot.check(br.getZ(), a1.getZ())) {
						forceVal += pot.calcF(intraMolecVec.length());
						potVal += pot.calcU(intraMolecVec.length());
					}
			}
			System.out.println((startVal + i*step) + "\t" + potVal + "\t" + forceVal);
		}
		double brBrDistAfter = JVector.distance(atoms[1].getPosition(), atoms[2].getPosition());
		double cBrDistAfter = JVector.distance(atoms[0].getPosition(), atoms[1].getPosition());
		double brBrForceAfter = ljp.calcF(brBrDistAfter);
		double cBrForceAfter = hooke.calcF(cBrDistAfter);
		System.out.println("Br-Br intermolecular distance/force before relaxing: " + 
				brBrDistBefore + "/" + brBrForceBefore);
		System.out.println("Br-Br intermolecular distance/force after relaxing: " + 
				brBrDistAfter + "/" + brBrForceAfter);
		System.out.println("C-Br intermolecular distance/force before relaxing: " + 
				cBrDistBefore + "/" + cBrForceBefore);
		System.out.println("C-Br intermolecular distance/force after relaxing: " + 
				cBrDistAfter + "/" + cBrForceAfter);
		
		System.out.println("New C-Br bond distance: " + newRestingPosition);
		System.out.println("Total Br-Br force: " + brBrforce.length());
		
	}
	
	public static void main(String[] args) {
		new PotentialTest().run();
	}
}
