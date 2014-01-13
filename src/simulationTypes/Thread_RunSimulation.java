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
import indexing.AlphabeticIndexingSystem;
import io.MyFile;
import io.MyObjectOutputStream;
import io.MyPrintStream;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Random;
import java.util.Vector;

import simulationTools.LatticeTools;
import defaultPackage.IdealTetrahedron;
import defaultPackage.Lattice;
import defaultPackage.LennardJonesLookup;
import defaultPackage.LennardJonesPotential;
import defaultPackage.Simulate;


public class Thread_RunSimulation implements Runnable {

	
	private File inputFolder;
	private int numThreadsToUse = 1;
	private SimulationToRun simulType;
	private Random r;
	private long randomNumberSeed = System.currentTimeMillis();
	
	/* SIMULATION OUTPUT VARIABLES */
	private File outputFolder, objOutputFolder, xyzOutputFolder, analysisOutputFolder, diffractionOutputFolder; 
	private volatile MyPrintStream mps_log;
	private volatile MyObjectOutputStream moos;
	private String initials = "EDD";
	private int notebookNumber = 4;
	private int pageNumber = 45;
	private AlphabeticIndexingSystem ais = new AlphabeticIndexingSystem("a");
	
	/* ----- */
	/* INPUT */
	private File[] initialObjFiles;
	private File[] modifiedObjFiles;
	private File[] xyzFiles;
	/* END INPUT */
	/* --------- */
	
	/* -------------------------- */
	/* INITIAL LATTICE PARAMETERS */
	private InitialLatticeTypes initLattice;
	private int numLatticesToMake = 1;
	private double a = 8.82;
	private int numUnitCellsPerAxis = 5;
	private String initLatticeObjFileExtension = ".init.obj";
	/* END INITIAL LATTICE PARAMETERS */
	/* ------------------------------ */
	
	/* ------------------------- */
	/* MODIFY LATTICE PARAMETERS */
	private ModifyLatticeTypes modLattice;
	private int startModifyIdx = 0;
	private int finishModifyIdx = 10;
	private String modLatticeObjFileExtension = ".mod.obj";
	private double minimumIntermolecularDistance = 0;
	
	/* interaction parameters */
	private LennardJonesLookup ljl;
	private double lookupPrecision = 0.001;
	
	private LennardJonesPotential ljp;
	private int Z1 = 35;
	private int Z2 = 35;
	private double Z1Z2_distance = 3.8;
	private double zeroPotentialDistance;
	private double wellDepth = 0.17;
	private String energyUnits = "eV";
	
	/* random walk parameters */
	private int numRandomWalksPerMolecule = 5;
	private int numStepsPerWalk = 256;
	private double maximumTranslationalMotionPerMove = 0.25;
	private double maximumRotationalAnglePerRotation = 1;
	private double energyCutoff = 0;
	private boolean useEnergyCutoff = false;
	private boolean disorderBetweenCycles = false;
	
	/* monte carlo orientational changes parameters */
	private int numTestsPerMolecule = 5;
	private boolean maxNumChangesEnabled = false;
	private int maxNumOrientationalChanges = 1000;
	/* END MODIFY LATTICE PARAMETERS */
	/* ----------------------------- */	
	
	/* ------------------------------ */
	/* CALCULATE DIFFRACTION PARAMETERS */

	/* END CALCULATE DIFFRACTION PARAMETERS */
	/* ------------------------------------ */
	
	/* -------------------------- */	
	/* ANALYZE_LATTICE PARAMETERS */

	/* END ANALYZE_LATTICE PARAMETERS */
	/* ------------------------------ */
	
	
	public Thread_RunSimulation(File outputRoot, SimulationToRun simulType) {
		makeOutputFolder(outputRoot);
		this.simulType = simulType;
		initStreams();
		initInputFileArrays();
	}
	public Thread_RunSimulation(File outputRoot, File inputFolder, SimulationToRun simulType) {
		makeOutputFolder(outputRoot);
		this.inputFolder = inputFolder;
		this.simulType = simulType;
		initStreams();
		initInputFileArrays();
	}
	private void initInputFileArrays() {
		if(inputFolder == null)
			return;
		
		File objInput = new File(inputFolder + File.separator + "obj");
		File xyzInput = new File(inputFolder + File.separator + "xyz"); 
		try {
			File[] files = objInput.listFiles();
			Vector<File> initFiles = new Vector<File>();
			Vector<File> modFiles = new Vector<File>();
			for(int i = 0; i < files.length; i++) {
				if(files[i].getName().contains(initLatticeObjFileExtension)) {
					initFiles.add(files[i]);
				} else if(files[i].getName().contains(modLatticeObjFileExtension)) {
					modFiles.add(files[i]);
				}
			}
			if(initFiles.size() > 0) {
				initialObjFiles = new File[initFiles.size()];
				initialObjFiles = initFiles.toArray(initialObjFiles);
			}
			if(modFiles.size() > 0) {
				modifiedObjFiles = new File[modFiles.size()];
				modifiedObjFiles = modFiles.toArray(modifiedObjFiles);
			}
		} catch(Exception e) {
			writeToLog("InitialLattice object file folder: " + objInput.getAbsolutePath() + " does not exist.");
			StackTraceElement[] ste = e.getStackTrace();
			writeToLog("Error: " + ste[0].toString());
			for(int i = 1; i < ste.length; i ++) {
				writeToLog("\tError: " + ste[i].toString());
			}
		}
		try {
			xyzFiles = new File(inputFolder + File.separator + "xyz").listFiles();
		} catch(Exception e) {
			writeToLog("xyz file folder: " + xyzInput.getAbsolutePath() + " does not exist.");
			StackTraceElement[] ste = e.getStackTrace();
			writeToLog("Error: " + ste[0].toString());
			for(int i = 1; i < ste.length; i ++) {
				writeToLog("\tError: " + ste[i].toString());
			}
		}
	}
	private void makeOutputFolder(File outputRoot) {
		outputFolder = new File(outputRoot + File.separator + getExptName());
		if(!outputFolder.mkdirs()) {
			ais.update();
			makeOutputFolder(outputRoot);
		}
		
	}
	private String getExptName() {
		return initials + "_" + notebookNumber + "-" + pageNumber + ais.getName();
	}
	private void initStreams() {
		outputFolder.mkdirs();
		mps_log = new MyPrintStream(new File(outputFolder + File.separator + getExptName() + ".log"));
		writeToLog("Working Directory: " + mps_log.getFile().getAbsolutePath());
		writeToLog("Simulation type: " + simulType.name());
	}
	/**
	 * 
	 * @param numToMake
	 * @param latticeType
	 */
	private void makeInitialLattice() {
		assert initLattice != null;
		writeToLog("Making " + numLatticesToMake + " initial lattices.");
		writeToLog("\tLattice Type: " + initLattice.name());
		writeToLog("\tLattice Parameters: ");
		writeToLog("\t\ta = " + a);
		writeToLog("\t\tnumUnitCellsPerAxis = " + numUnitCellsPerAxis);
		switch(initLattice) {
		case RANDOM_SIXFOLD:
			Lattice l = new Lattice(numUnitCellsPerAxis, a);
			l.setRandom(r);
			IdealTetrahedron[][][] lattice;
			InitialLattice init;
			MyObjectOutputStream moos;
			MyPrintStream mps;
			writeToLog("\tcalling l.makeTetragonalFCCLattice() " + numLatticesToMake + " times.");
			
			for(int i = 0; i < numLatticesToMake; i++) {
				writeToLog("\tLattice " + (i+1) + " of " + numLatticesToMake);
				
				lattice = l.makeTetragonalFCCLattice();
				writeToLog("\t\t made tetragonal (six options) FCC Lattice");
				
				init = new InitialLattice();
				init.setA(a);
				init.setLattice(lattice);
				init.setNumUnitCellsPerAxis(numUnitCellsPerAxis);
				init.setInitLatticeObjFileExtension(initLatticeObjFileExtension);
				
				moos = new MyObjectOutputStream(new File(objOutputFolder + File.separator + i + ".init.obj"));
				moos.writeObject(init);
				moos.close();
				writeToLog("\t\tObject written to file.");
				
				mps = new MyPrintStream(new File(xyzOutputFolder + File.separator + i + ".init.xyz"));
				LatticeTools.printLatticeXYZ(lattice, mps.getPrintStream(), a);
				mps.close();
				writeToLog("\t\tLattice written to file.");
				mps_log.flush();
			}
			break;
		}
		writeToLog("Making " + numLatticesToMake + " initial lattices complete.");
	}
	private void modifyLattice() {
		assert initLattice != null;
		writeToLog("Modifying " + (finishModifyIdx - startModifyIdx) + " initial lattices.");
		writeToLog("\tModification Type: " + modLattice.name());
		writeToLog("\tLattice Parameters: ");
		writeToLog("\t\ta = " + a);
		writeToLog("\t\tnumUnitCellsPerAxis = " + numUnitCellsPerAxis);
		
		zeroPotentialDistance = Z1Z2_distance / Math.pow(2., (1./6.));
		ljp = new LennardJonesPotential(Z1, Z2, zeroPotentialDistance, wellDepth);
		ljl = new LennardJonesLookup(lookupPrecision, ljp);
		writeToLog("LennardJonesPotential:");
		writeToLog("\tZ1: " + Z1);
		writeToLog("\tZ2: " + Z2);
		writeToLog("\tWell depth (epsilon): " + wellDepth);
		writeToLog("\tZ1-Z2 distance: " + Z1Z2_distance);
		writeToLog("\tZero potential distance (sigma): " + zeroPotentialDistance);
		writeToLog("LennardJonesLookup precision: " + lookupPrecision);
		
		switch(modLattice) {
		case RANDOM_WALK:
			InitialLattice init = null;
			ModifiedLattice mod;
			MyPrintStream mps;
			writeToLog("\tPerforming Random Walk energy minimization.");
			writeToLog("\t\tnumRandomWalksPerMolecule: " + numRandomWalksPerMolecule);
			writeToLog("\t\tnumStepsPerWalk: " + numStepsPerWalk);
			writeToLog("\t\tmaximumTranslationalMotionPerMove: " + maximumTranslationalMotionPerMove);
			writeToLog("\t\tmaximumRotationalAnglePerRotation: " + maximumRotationalAnglePerRotation);
			writeToLog("\t\tUsing an energy cutoff value (" + useEnergyCutoff + ") of: " + energyCutoff);
			writeToLog("\t\tDisordering between cycles: " + disorderBetweenCycles);
			
			for(int i = startModifyIdx; i < finishModifyIdx; i++) {
				writeToLog("\tModifying lattice " + (i) + ". Started at: " + startModifyIdx + " and ending at: " + finishModifyIdx);
				try {
					init = LatticeTools.readInInitialLattice(initialObjFiles[i]);
				} catch(IOException e) {
					e.printStackTrace();
				} catch(ClassNotFoundException e) {
					e.printStackTrace();
				}
				writeToLog("\t\tSuccessfully read in InitialLattice object");
				Simulate simul = new Simulate(ljl, minimumIntermolecularDistance, init.getA());
				simul.setMaxTorque(maximumRotationalAnglePerRotation);
				simul.setMaxTrans(maximumTranslationalMotionPerMove);
				simul.setLattice3(init.getLattice());
				simul.setRandom(r);
				writeToLog("\tStarting random walk.");
				simul.calculateTotalEnergy();
				writeToLog("\t\tEnergy before random walk starts: " + simul.getTotalEnergy() + " " + energyUnits);
				for(int j = 0; j < numRandomWalksPerMolecule; j++) {
					simul.walk(numStepsPerWalk, 1, energyCutoff, useEnergyCutoff, disorderBetweenCycles);
					simul.calculateTotalEnergy();
					writeToLog("\t\tEnergy after walk " + (j+1) + " of " + numRandomWalksPerMolecule + ": " + simul.getTotalEnergy() + " " + energyUnits);
				}
				mod = new ModifiedLattice(init);
				setModLatticeParams(mod);
				
				moos = new MyObjectOutputStream(new File(objOutputFolder + File.separator + i + ".mod.obj"));
				moos.writeObject(mod);
				moos.close();
				writeToLog("\t\tModified lattice object written to file.");
				
				mps = new MyPrintStream(new File(xyzOutputFolder + File.separator + i + ".mod.xyz"));
				LatticeTools.printLatticeXYZ(mod.getModifiedLattice(), mps.getPrintStream(), a);
				mps.close();
				writeToLog("\t\tLattice written to file.");
				mps_log.flush();
			}
			break;
		case MONTE_CARLO_ORIENTATIONAL:
			
			break;
		}
		writeToLog("Making " + numLatticesToMake + " initial lattices complete.");
	}
	private void setModLatticeParams(ModifiedLattice mod) {
		mod.setModLattice(modLattice);
		mod.setStartModifyIdx(startModifyIdx);
		mod.setFinishModifyIdx(finishModifyIdx);
		mod.setModifiedLatticeObjFileName(modLatticeObjFileExtension);
		mod.setMinimumIntermolecularDistance(minimumIntermolecularDistance);
		
		/* interaction parameters */
		mod.setLookupPrecision(lookupPrecision);
		mod.setLjp(ljp);
		mod.setZ1(Z1);
		mod.setZ2(Z2);
		mod.setZ1Z2_distance(Z1Z2_distance);
		mod.setZeroPotentialDistance(zeroPotentialDistance);
		mod.setWellDepth(wellDepth);
		mod.setEnergyUnits(energyUnits);
		
		/* random walk parameters */
		mod.setNumRandomWalksPerMolecule(numRandomWalksPerMolecule);
		mod.setNumStepsPerWalk(numStepsPerWalk);
		mod.setMaximumTranslationalMotionPerMove(maximumTranslationalMotionPerMove);
		mod.setMaximumRotationalAnglePerRotation(maximumRotationalAnglePerRotation);
		mod.setEnergyCutoff(energyCutoff);
		mod.setUseEnergyCutoff(useEnergyCutoff);
		mod.setDisorderBetweenCycles(disorderBetweenCycles);

		/* monte carlo orientational changes parameters */
		mod.setNumTestsPerMolecule(numTestsPerMolecule);
		mod.setMaxNumChangesEnabled(maxNumChangesEnabled);
		mod.setMaxNumOrientationalChanges(maxNumOrientationalChanges);
		
		mod.setRandomNumberSeed(randomNumberSeed);
		
		
	}
	private void calcDiffraction() {
		
	}
	
	private void analyzeLattice() {
		
	}

	private void initObjOut() {
		objOutputFolder = new File(outputFolder.getAbsoluteFile() + File.separator + "obj");
		objOutputFolder.mkdirs();
	}
	private void initXYZOut() {
		xyzOutputFolder = new File(outputFolder.getAbsoluteFile() + File.separator + "xyz");
		xyzOutputFolder.mkdirs();
	}
	public void run() {
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		writeToLog("RunSimulation.run() started at: " + dateFormat.format(new Date()));
		r = new Random(randomNumberSeed);
		writeToLog("Random number seed: " + randomNumberSeed);
		writeToLog("Using: " + numThreadsToUse + " threads.");
		switch(simulType) {
		case INITIAL_LATTICE:
			initObjOut();
			initXYZOut();
			makeInitialLattice();
			break;
		case MODIFY_LATTICE:
			initObjOut();
			initXYZOut();
			modifyLattice();
			break;
		case CALC_DIFFRACTION:
			calcDiffraction();
			break;
		case ANALYZE_LATTICE:
			analyzeLattice();
			break;
		}
		writeToLog("RunSimulation.run() complete at: " + dateFormat.format(new Date()));
	}
	enum SimulationToRun {
		INITIAL_LATTICE,
		MODIFY_LATTICE,
		CALC_DIFFRACTION,
		ANALYZE_LATTICE,
		;
	}
	
	private synchronized void writeToLog(String toWrite) {
		System.out.println(toWrite);
		mps_log.println(toWrite);
	}
	public static void main(String[] args) {
		String root = "C:\\$tempResearch\\";
		String user = "Dill";
		int notebookNumber = 4;
		int pageNumber = 83;
		
		File outputFile = new File(root + user + "-" + notebookNumber + File.separator + pageNumber);
		File inputFile = new File("C:\\$tempResearch\\Dill-4\\82\\EDD_4-45p");
		
		Thread_RunSimulation simul = new Thread_RunSimulation(outputFile, inputFile, SimulationToRun.INITIAL_LATTICE);
		simul.numUnitCellsPerAxis = 10;
		simul.numLatticesToMake = 100;
		simul.initLattice = InitialLatticeTypes.RANDOM_SIXFOLD;
		simul.modLattice = ModifyLatticeTypes.RANDOM_WALK;
		simul.run();
	}
}
