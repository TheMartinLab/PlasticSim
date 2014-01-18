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
package simulationTypes;

import io.MyPrintStream;

import java.io.File;

import simulationTools.LatticeTools;
import chemistry.JAtomTools;
import defaultPackage.IdealTetrahedron;
import defaultPackage.JAtom;
import defaultPackage.JVector;
import defaultPackage.Molecule;
import defaultPackage.Potential;
import defaultPackage.PotentialLookup;
import defaultPackage.Simulate;
import formatting.FormatDouble;

public class Thread_DistortTetrahedra extends Thread {

	private DistortTetrahedra distort;
	private Simulate simul;
	private int checkEnergyEvery;
	private FormatDouble fd;
	private int timeStep;
	public Thread_DistortTetrahedra(DistortTetrahedra distort) {
		this.distort = distort;
		fd = new FormatDouble(FormatDouble.DecimalPlaces.THREE);
	}

	public void close() {
		distort.close();
	}

	@Override
	public void run() {
		/* LOG PARAMETERS */
		distort.logParameters();
		distort.getLog().println("\n\n-- RUNTIME INFORMATION --\n\n");
		/* LATTICE PARAMETERS */
		IdealTetrahedron[][][] lattice = distort.getLattice().getLattice();
		simul = new Simulate(null, 0, distort.getLattice().getA());
		simul.setLattice3(lattice);
		for(PotentialLookup potLook : distort.getInterMolecPotentials()) {
			simul.addPotentialLookup(potLook);
		}
		Molecule[] mols = simul.to1dArray();
		Molecule[][] surrounding = simul.getSurrounding();
		log("Surrounding molecule matrix constructed.", "\t");
		
		/* TIME PARAMETERS */
		double time = 0;
		
		/* ENERGY PARAMETERS */
		double prevLatticeEnergy = 0;
		double curLatticeEnergy = 0;
		double prevMaxMoved = 0;
		double deltaE = 0;
		boolean hasCheckedEnergy = false;
		boolean continueLoop = true;
		
		/* INITIAL SIMULATION PARAMETERS */
		simul.calculateTotalEnergy();
		curLatticeEnergy = simul.getTotalEnergy();
		JVector[][] forces = new JVector[mols.length][5];
		log("Initial lattice energy\t" + curLatticeEnergy, "\t");
		timeStep = 0;
		
		/* SET UP MOVIE OUTPUT */
		MyPrintStream mps_movie = null;
		if(distort.isOutputtingMovie()) { 
			File mps_file = new File(distort.getOutputFolder() + File.separator + "movies");
			mps_file.mkdirs();
			mps_file = new File(mps_file + File.separator + distort.getOutputFilePrefix() + "movie.xyz");
			mps_movie = new MyPrintStream(mps_file); 
		}
		do {
			timeStep++;
			log("Time step # " + timeStep + " = " + Math.rint(time/distort.getTimeStep()) * distort.getTimeStep(), "");
			/* DO ENERGY MINIMIZATION STUFF */
			zeroForces(forces);
			calcForcesOnMolecules(mols, surrounding, forces);
			prevMaxMoved = move(mols, forces, prevMaxMoved);
			time += distort.getTimeStep();

			/* OUTPUT XYZs */
			if(distort.isOutputtingXYZs() && 
					timeStep % distort.getOutputXYZsEveryThisManyTimeSteps() == 0) {
				log("Outputting XYZ", "\t");
				outputXYZ(mols, time);
			}
			
			/* OUTPUT MOVIE */
			if(distort.isOutputtingMovie() && 
					timeStep % distort.getOutputMovieXYZsEveryThisManyTimeSteps() == 0) {
				log("Outputting Movie XYZ", "\t");
				outputMovie(mols, mps_movie);
			}
			
			/* CHECK LATTICE ENERGY */
			if(timeStep % checkEnergyEvery == 0) {
				prevLatticeEnergy = curLatticeEnergy;
				simul.calculateTotalEnergy();
				curLatticeEnergy = simul.getTotalEnergy();
				log("Lattice energy at " + time + ":\t" + curLatticeEnergy, "\t");
				deltaE = curLatticeEnergy - prevLatticeEnergy;
				hasCheckedEnergy = true;
			}
			continueLoop = !endConditionsAreMet(deltaE, prevMaxMoved);
//			else
//				continueLoop = !endConditionsAreMet(deltaE, prevMaxMoved);
			
		} while(continueLoop);
		
		/* SHUT DOWN MOVIE PRINT STREAM */
		if(distort.isOutputtingMovie()) {
			mps_movie.flush();
			mps_movie.close();
		}
	}
	private void log(String toLog, String prefix) {
		System.out.println(prefix + distort.getOutputFilePrefix() + ": " + toLog);
		distort.getLog().println(prefix + toLog);
	}
	private boolean endConditionsAreMet(double deltaE, double prevMaxMoved) {
		boolean end = false;
		if(!end && distort.isEndingBecauseOfEnergy()) {
			end = Math.abs(deltaE) < distort.getEndDeltaE();
			if(end)	
				log("END CONDITION: ENDING BECAUSE MINIMUM DELTA E HAS BEEN MET: " + deltaE + "<" + distort.getEndDeltaE(), "");
		}
		if(!end && distort.isEndingBecauseOfLargeMovement()) {
			end = prevMaxMoved > distort.getLargeEndMovement();
			if(end)	
				log("END CONDITION: ENDING BECAUSE MOVEMENT EXCEEDS MAXIMUM THRESHOLD: " + prevMaxMoved + "<" + distort.getLargeEndMovement(), "");
		}
		if(!end && distort.isEndingBecauseOfSmallMovement()) {
			end = prevMaxMoved < distort.getSmallEndMovement();
			if(end)	
				log("END CONDITION: ENDING BECAUSE MOVEMENT IS LESS BELOW MINIMUM THRESHOLD: " + prevMaxMoved + "<" + distort.getSmallEndMovement(), "");
		}
		if(!end && distort.isEndingBecauseOfTimeSteps()) {
			end = timeStep > distort.getEndTimeStep();
			if(end)	
				log("END CONDITION: ENDING BECAUSE > " + distort.getEndTimeStep() + " TIME STEPS HAVE OCCURED", "");
		}
		
		return end;
	}
	private JAtom[] getAtoms(Molecule[] mols) {
		JAtom[] atoms = new JAtom[mols.length * 5];
		int atomsPerMolecule = 5;
		for(int i = 0; i < mols.length; i++) {
			JAtom[] perMol = mols[i].getAtoms();
			for(int a = 0; a < perMol.length; a++) {
				atoms[i*atomsPerMolecule + a] = perMol[a];
			}
		}
		return atoms;
	}
	private void outputMovie(Molecule[] mols, MyPrintStream mps_movie) {
		mps_movie.println(mols.length  * 5);
		mps_movie.println();
		JAtom[] atoms = getAtoms(mols);
		String line = "";
		for(int i = 0; i < atoms.length; i++) {
			line += atoms[i].getZ() + "\t" + atoms[i].getPosition().toTabString() + "\n";
			/* ONLY PRINT EVERY 1000 LINES TO INCREASE OUTPUT SPEED */
			if(i % 1000 == 0) {
				mps_movie.print(line);
				line = "";
			}
		}
		mps_movie.print(line);
		mps_movie.flush();
	}
	private void outputXYZ(Molecule[] mols, double timeStep) {
		MyPrintStream mps = new MyPrintStream(new File(distort.getOutputFolder() + File.separator + distort.getOutputFilePrefix() + "." + fd.format(timeStep) + ".xyz"));
		JAtom[] atoms = getAtoms(mols);
		LatticeTools.printLatticeXYZ(atoms, mps.getPrintStream());
		mps.flush();
		mps.close();
	}
	private double move(Molecule[] mols, JVector[][] forces, double prevMaxMoved) {
		JVector force;
		double mass;
		double forceMagnitude;
		double timeStep = distort.getTimeStep();
		JVector maxMovement = new JVector();
		for(int i = 0; i < mols.length; i++) {
			JAtom[] atoms = mols[i].getAtoms();
			for(int a = 0; a < atoms.length; a++) {
				force = forces[i][a];
				mass = JAtomTools.getMass(atoms[a].getZ());
				forceMagnitude = 0.5 * mass * timeStep * timeStep;
				force = JVector.multiply(force, forceMagnitude);
				if(!distort.isEndingBecauseOfLargeMovement() && force.length() > distort.getLargeEndMovement())
					force = JVector.multiply(force.unit(), distort.getLargeEndMovement());
				if(maxMovement.length() < force.length())
					maxMovement = force;
				atoms[a].getPosition().add_inPlace(force);
			}
		}
		double maxAmountMoved = maxMovement.length();
		
		log("Max amount moved\t" + maxAmountMoved, "\t");
		
//		if(prevMaxMoved != 0 && maxAmountMoved > prevMaxMoved) {
//			System.out.println("\tPrevious time step: " + distort.getTimeStep());
//			distort.setTimeStep(distort.getTimeStep() * 0.9);
//			System.out.println("\tChanging time step to: " + distort.getTimeStep());
//		} 
		return maxAmountMoved;
	}
	private void zeroForces(JVector[][] forces) {
		for(int i = 0; i < forces.length; i++)
			for(int j = 0; j < forces[i].length; j++)
				if(forces[i][j] == null)
					forces[i][j] = new JVector();
				else
					forces[i][j].zero();
	}
	private void calcForcesOnMolecules(Molecule[] mols, Molecule[][] surroundings, JVector[][] forces) {
		int Z1, Z2;
		JVector forceVec = null, curPos = null, surrPos = null, translationVec = null, curCenter = null, surrCenter = null;
		double force;
		boolean haveBeenMoved = false;
		
		Molecule curMol = null, surrMol = null;
		JAtom curArr[], surrArr[], curAtom, surrAtom;
		for(int i = 0; i < mols.length; i++) {
			curMol = mols[i];
			curArr = curMol.getAtoms();
			for(int a = 0; a < curArr.length; a++) {
				curAtom = curArr[a];
				Z1 = curAtom.getZ();
				curPos = curAtom.getPosition();
				if(Z1 == 35) {
					/* ************************** */
					/* INTERMOLECULAR CALCULATION */
					/* ************************** */
					for(int j = 0; j < surroundings[i].length; j++) {
						surrMol = surroundings[i][j];
						translationVec = areTooFarApart(curMol.getCenter(), surrMol.getCenter());
						
						if(translationVec != null) {
							haveBeenMoved = true;
							surrMol.translate(translationVec);
						}
						surrArr = surrMol.getAtoms();
						for(int b = 0; b < surrArr.length; b++) {
							surrAtom = surrArr[b];
							Z2 = surrAtom.getZ();
							surrPos = surrAtom.getPosition();
							if(Z2 == 35)
								for(PotentialLookup potLook : distort.getInterMolecPotentials()) {
									if(potLook.check(Z1, Z2)) {
										forceVec = JVector.subtract(curPos, surrPos);
										force = potLook.lookupForce(forceVec.length());
										forceVec = JVector.multiply(forceVec.unit(), force/2.);
										forces[i][a].add_inPlace(forceVec);
									}
								}
						}
						if(haveBeenMoved) {
							haveBeenMoved = false;
							translationVec.multiply(-1);
							surrMol.translate(translationVec);
							translationVec = null;
						}
					}
				}
				/* ************************** */
				/* INTRAMOLECULAR CALCULATION */
				/* ************************** */
				for(int b = a; b < curArr.length; b++) {
					surrAtom = curArr[b];
					if(curAtom == surrAtom)
						continue;
					Z2 = surrAtom.getZ();
					surrPos = surrAtom.getPosition();
					for(PotentialLookup potLook : distort.getIntraMolecPotentials()) {
						if(potLook.check(Z1, Z2)) {
							forceVec = JVector.subtract(curPos, surrPos);
							force = potLook.lookupForce(forceVec.length());
							forceVec = JVector.multiply(forceVec.unit(), force/2.);
							forces[i][a].add_inPlace(forceVec);
						}
					}
				}
			}
		}
	}
	private JVector areTooFarApart(JAtom center, JAtom surrounding) {
		boolean areTooFarApart = false;
		double maxDist = distort.getIntermolecularMaxDist();
		JVector translate = new JVector();
		int numUnitsPerAxis = distort.getLattice().getNumUnitCellsPerAxis();
		double latticeParam = distort.getLattice().getA();
		double positiveTranslate = numUnitsPerAxis * latticeParam;
		double negativeTranslate = -1 * positiveTranslate;
		
		double x = center.getPosition().i - surrounding.getPosition().i;
		double y = center.getPosition().j - surrounding.getPosition().j;
		double z = center.getPosition().k - surrounding.getPosition().k;
		if(Math.abs(x) > maxDist) {
			areTooFarApart = true;
			if(x > 0)
				translate.i = positiveTranslate;
			else
				translate.i = negativeTranslate;
		}
		if(Math.abs(y) > maxDist) {
			areTooFarApart = true;
			if(y > 0)
				translate.j = positiveTranslate;
			else
				translate.j = negativeTranslate;
		}
		if(Math.abs(z) > maxDist) {
			areTooFarApart = true;
			if(z > 0)
				translate.k = positiveTranslate;
			else
				translate.k = negativeTranslate;
		}
		if(!areTooFarApart)
			translate = null;
		
		return translate;
	}
	public DistortTetrahedra getDistort() { return distort; }
	public void setDistort(DistortTetrahedra distort) { this.distort = distort; }
	public Simulate getSimul() { return simul; }
	public void setSimul(Simulate simul) { this.simul = simul; }
	public int getCheckEnergyEvery() { return checkEnergyEvery; }
	public void setCheckEnergyEvery(int checkEnergyEvery) { this.checkEnergyEvery = checkEnergyEvery; }
}
