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

import io.MyObjectOutputStream;
import io.MyPrintStream;

import java.io.File;
import java.io.IOException;
import java.util.Observable;
import java.util.Random;

import simulationTools.LatticeTools;
import defaultPackage.LennardJonesPotential;
import defaultPackage.PotentialLookup;
import defaultPackage.Simulate;

public class Thread_ModifyLattice extends Observable implements Runnable {
	private ModifiedLattice mod;
	private String output;
	public final static String COMPLETE = "Modify Lattice Complete.";
	private MyPrintStream mps_log;
	private Thread t;
	
	public Thread_ModifyLattice(ModifiedLattice mod, MyPrintStream mps_log) {
		this.mod = mod;
		this.mps_log = mps_log;
		output = "\nThread " + mod.getThreadIdx() + " output start.";
		writeToLog("\tThread " + mod.getThreadIdx() + " output start.");
		t = new Thread(this);
//		t.
	}
	
	private void writeToLog(String s) {
		mps_log.println("Thread: " + mod.getThreadIdx() + " " + s);
	}
	public void run() {
		writeToLog("\tModification Type: " + mod.getModLattice().name());
		writeToLog("\tLattice Parameters: ");
		writeToLog("\t\ta = " + mod.getA());
		writeToLog("\t\tnumUnitCellsPerAxis = " + mod.getNumUnitCellsPerAxis());
		
		double zeroPotentialDistance = mod.getZ1Z2_distance() / Math.pow(2., (1./6.));
		mod.setZeroPotentialDistance(zeroPotentialDistance);
		LennardJonesPotential ljp = new LennardJonesPotential(mod.getZ1(), mod.getZ2(), mod.getZeroPotentialDistance(), mod.getWellDepth());
		PotentialLookup ljl = new PotentialLookup(mod.getLookupPrecision(), ljp);
		writeToLog("LennardJonesPotential:");
		writeToLog("\tZ1: " + mod.getZ1());
		writeToLog("\tZ2: " + mod.getZ2());
		writeToLog("\tWell depth (epsilon): " + mod.getWellDepth());
		writeToLog("\tZ1-Z2 distance: " + mod.getZ1Z2_distance());
		writeToLog("\tZero potential distance (sigma): " + mod.getZeroPotentialDistance());
		writeToLog("LennardJonesLookup precision: " + mod.getLookupPrecision());

		Lattices init = null;
		MyPrintStream mps;
		Simulate simul;
		MyObjectOutputStream moos;
		switch(mod.getModLattice()) {
		case RANDOM_WALK:
			System.out.println("Random walk energy minimization algorithm:");
			writeToLog("\tPerforming Random Walk energy minimization.");
			writeToLog("\t\tnumRandomWalksPerMolecule: " + mod.getNumRandomWalksPerMolecule());
			writeToLog("\t\tnumStepsPerWalk: " + mod.getNumStepsPerWalk());
			writeToLog("\t\tmaximumTranslationalMotionPerMove: " + mod.getMaximumTranslationalMotionPerMove());
			writeToLog("\t\tmaximumRotationalAnglePerRotation: " + mod.getMaximumRotationalAnglePerRotation());
			writeToLog("\t\tUsing an energy cutoff value (" + mod.isUseEnergyCutoff() + ") of: " + mod.getEnergyCutoff());
			writeToLog("\t\tDisordering between cycles: " + mod.isDisorderBetweenCycles());
			
			try {
				init = LatticeTools.readInLattice(new File(mod.getInputLatticeObjFileName()));
			} catch(IOException e) {
				e.printStackTrace();
			} catch(ClassNotFoundException e) {
				e.printStackTrace();
			}
			mod.setCommonLatticeParams(init);
			writeToLog("\t\tSuccessfully read in InitialLattice object from: " + mod.getInputLatticeObjFileName());
			simul = new Simulate(ljl, mod.getMinimumIntermolecularDistance(), mod.getA());
			simul.setMaxTorque(mod.getMaximumRotationalAnglePerRotation());
			simul.setMaxTrans(mod.getMaximumTranslationalMotionPerMove());
			simul.setLattice3(mod.getLattice());
			simul.setRandom(new Random(mod.getRandomNumberSeed()));
			writeToLog("\tStarting random walk.");
			simul.calculateTotalEnergy();
			writeToLog("\t\t" + mod.getThreadIdx() + ": Energy before random walk starts: " + simul.getTotalEnergy() + " " + mod.getEnergyUnits());
			for(int j = 0; j < mod.getNumRandomWalksPerMolecule(); j++) {
				simul.walk(mod.getNumStepsPerWalk(), 1, mod.getEnergyCutoff(), mod.isUseEnergyCutoff(), mod.isDisorderBetweenCycles());
				simul.calculateTotalEnergy();
				writeToLog("\t\t" + mod.getThreadIdx() + ": Energy after walk " + (j+1) + " of " + mod.getNumRandomWalksPerMolecule() + ": " + simul.getTotalEnergy() + " " + mod.getEnergyUnits());
			}
			mod.setLattice(init.getLattice());
			moos = new MyObjectOutputStream(new File(mod.getOutputLatticeObjFileName()));
			moos.writeObject(mod);
			moos.close();
			writeToLog("\t\tModified lattice object written to: " + mod.getOutputLatticeObjFileName());
			
			mps = new MyPrintStream(new File(mod.getOutputLatticeXYZFileName()));
			LatticeTools.printLatticeXYZ(mod.getLattice(), mps.getPrintStream(), mod.getA());
			mps.close();
			writeToLog("\t\tModified Lattice written to: " + mod.getOutputLatticeXYZFileName());
			setChanged();
			output += "\nThread " + mod.getThreadIdx() + " output stop.";
			notifyObservers("Thread: " + mod.getThreadIdx() + " complete.\nOutput:\n" + output);
			break;
		case MONTE_CARLO_ORIENTATIONAL:
			System.out.println("Monte carlo orientational change energy minimization algorithm:");
			writeToLog("\tPerforming Monte Carlo Orientational energy minimization.");
			writeToLog("\t\tnumTestsPerMolecule: " + mod.getNumTestsPerMolecule());
			writeToLog("\t\tmaxNumChangesEnabled: " + mod.isMaxNumChangesEnabled());
			writeToLog("\t\tmaxNumOrientationalChanges: " + mod.getMaxNumOrientationalChanges());
			
			try {
				init = LatticeTools.readInLattice(new File(mod.getInputLatticeObjFileName()));
			} catch(IOException e) {
				e.printStackTrace();
			} catch(ClassNotFoundException e) {
				e.printStackTrace();
			}
			mod.setCommonLatticeParams(init);
			writeToLog("\t\tSuccessfully read in InitialLattice object from: " + mod.getInputLatticeObjFileName());
			simul = new Simulate(ljl, mod.getMinimumIntermolecularDistance(), mod.getA());

			simul.setLattice3(mod.getLattice());
			simul.setRandom(new Random(mod.getRandomNumberSeed()));
			writeToLog("\tStarting Monte Carlo Orientational Changes.");
			simul.calculateTotalEnergy();

			writeToLog("\t\tEnergy beforemonte carlo orientational changes: " + simul.getTotalEnergy() + " " + mod.getEnergyUnits());
			System.out.println("\t\tEnergy beforemonte carlo orientational changes: " + simul.getTotalEnergy() + " " + mod.getEnergyUnits());
			
			int[] monteCarloInfo = simul.monteCarloOrientation(mod.getNumTestsPerMolecule());
			
			simul.calculateTotalEnergy();
			writeToLog("\t\tEnergy after monte carlo orientational changes: " + simul.getTotalEnergy() + " " + mod.getEnergyUnits());
			System.out.println("\t\tEnergy after monte carlo orientational changes: " + simul.getTotalEnergy() + " " + mod.getEnergyUnits());

			writeToLog("\t\tNumber of orientational tests: " + monteCarloInfo[0]);
			writeToLog("\t\tNumber of orientational changes: " + monteCarloInfo[1]);
			
			mod.setLattice(init.getLattice());
			moos = new MyObjectOutputStream(new File(mod.getOutputLatticeObjFileName()));
			moos.writeObject(mod);
			moos.close();
			writeToLog("\t\tModified lattice object written to: " + mod.getOutputLatticeObjFileName());
			
			mps = new MyPrintStream(new File(mod.getOutputLatticeXYZFileName()));
			LatticeTools.printLatticeXYZ(mod.getLattice(), mps.getPrintStream(), mod.getA());
			mps.close();
			writeToLog("\t\tModified Lattice written to: " + mod.getOutputLatticeXYZFileName());
			setChanged();
			output += "\nThread " + mod.getThreadIdx() + " output stop.";
			notifyObservers("Thread: " + mod.getThreadIdx() + " complete.\nOutput:\n" + output);
			break;
		}
		
		setChanged();
		notifyObservers(COMPLETE);
	}
}
