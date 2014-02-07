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
import io.MyObjectOutputStream;
import io.MyPrintStream;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Observable;
import java.util.Observer;
import java.util.Random;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import simulationTools.DiffractionTools;
import simulationTools.LatticeTools;
import simulationTools.WatchForNewFile;
import defaultPackage.IdealTetrahedron;
import defaultPackage.JPixel;
import defaultPackage.JVector;
import defaultPackage.Lattice;
import defaultPackage.LennardJonesPotential;
import defaultPackage.PotentialLookup;
import defaultPackage.Simulate;

public class RunSimulationThread implements Observer {

	/* ************ */
	/* CONSTRUCTORS */
	/* ************ */
	
	public RunSimulationThread(File outputRoot) {
		this.outputRoot = outputRoot;
	}
	public RunSimulationThread(File outputRoot, File inputFolder) {
		this.outputRoot = outputRoot;
		this.inputFolder = inputFolder;
	}
	
	/* *********************************** */
	/* INITIALIZATION AND INPUT PARAMETERS */
	/* *********************************** */
	private File[] objFiles;
	private File[] xyzFiles;
	private File inputFolder;
	private File outputRoot;
	private int numThreadsToUse = 1;
	private int numDiffractionThreadsToUse = 1; // this must be set to 1 for now. otherwise the calculation doesn't work. it runs, but it doesn't calculate correctly.
	private SimulationToRun simulType;
	private Random r;
	private long randomNumberSeed = System.currentTimeMillis();
	private ExecutorService es_simul, es_diffraction;
	private DateFormat df = DateFormat.getInstance();
	enum SimulationToRun {
		INITIAL_LATTICE,
		MODIFY_LATTICE,
		CALC_DIFFRACTION_OBJECT,
		CALC_DIFFRACTION_XYZ,
		ANALYZE_LATTICE,
		;
	}
	
	/* ******************************** */
	/* INITIALIZATION AND INPUT METHODS */
	/* ******************************** */
	private void initExecutors() {
		es_simul = Executors.newFixedThreadPool(numThreadsToUse);
		es_diffraction = Executors.newFixedThreadPool(numDiffractionThreadsToUse);
	}
	private void initInputFileArrays() {
		if(inputFolder == null)
			return;
		
		File objInput = new File(inputFolder + File.separator + "obj");
		File xyzInput = new File(inputFolder + File.separator + "xyz"); 
		try {
			objFiles = objInput.listFiles();
			
			if(finishModifyIdx == 0 || finishModifyIdx > objFiles.length)
				finishModifyIdx = objFiles.length;
				
		} catch(Exception e) {
			writeToLog("InitialLattice object file folder: " + objInput.getAbsolutePath() + " does not exist.");
			StackTraceElement[] ste = e.getStackTrace();
			try {
				writeToLog("Error: " + ste[0].toString());
				for(int i = 1; i < ste.length; i ++)
					writeToLog("\tError: " + ste[i].toString());
				}
			catch(ArrayIndexOutOfBoundsException ae) {
				// Do nothing...
			}
		}
		try {
			if(xyzFiles == null)
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
	private void initStreams() {
		if(outputFolder == null)
			makeOutputFolder(outputRoot);
			
		outputFolder.mkdirs();
		
		mps_log = new MyPrintStream(new File(outputFolder + File.separator + getExptName() + ".log"));
		
		writeToLog("Working Directory: " + mps_log.getFile().getAbsolutePath());
		writeToLog("Simulation type: " + simulType.name());
	}
	private void closeStreams() {
		mps_log.close();
	}
	private void makeOutputFolder(File outputRoot) {
		switch(simulType) {
		case ANALYZE_LATTICE:
		case CALC_DIFFRACTION_OBJECT:
			outputFolder = inputFolder;
			break;
		case INITIAL_LATTICE:
		case MODIFY_LATTICE:
			outputFolder = new File(outputRoot + File.separator + getExptName());
			if(!outputFolder.mkdirs()) {
				ais.update();
				makeOutputFolder(outputRoot);
			}
			break;
		default:
			break;
		}
	}
	private String getExptName() {
		return initials + "_" + notebookNumber + "-" + pageNumber + ais.getName();
	}
	private void initObjOut() {
		objOutputFolder = new File(outputFolder.getAbsoluteFile() + File.separator + "obj");
		objOutputFolder.mkdirs();
	}
	private void initXYZOut() {
		xyzOutputFolder = new File(outputFolder.getAbsoluteFile() + File.separator + "xyz");
		xyzOutputFolder.mkdirs();
	}
	/* ************************************ */
	/* END INITIALIZATION AND INPUT METHODS */
	/* ************************************ */
	
	
	/* ***************** */
	/* OUTPUT PARAMETERS */
	/* ***************** */
	private File outputFolder, objOutputFolder, xyzOutputFolder, analysisOutputFolder, diffractionOutputFolder;
	private File outputAllTo;
	private volatile MyPrintStream mps_log;
	private volatile MyObjectOutputStream moos;
	private String initials = "EDD";
	private int notebookNumber;
	private int pageNumber;
	private AlphabeticIndexingSystem ais = new AlphabeticIndexingSystem("a");
	
	/* ************** */
	/* OUTPUT METHODS */
	/* ************** */
	private synchronized void writeObject(Object obj, File fName) {
		moos = new MyObjectOutputStream(fName);
		moos.writeObject(obj);
		moos.close();
		writeToLog("\t\t" + obj.getClass() + " written to file.");
	}
	private synchronized void writeToLog(String toWrite) {
		System.out.println(df.format(new Date()) + "\t" + toWrite);
		mps_log.println(df.format(new Date()) + "\t" + toWrite);
	}
	
	/* ****************** */
	/* END OUTPUT METHODS */
	/* ****************** */
	/* ******************************* */
	/* MAKE INITIAL LATTICE PARAMETERS */
	/* ******************************* */
	private InitialLatticeTypes initLattice;
	private int numLatticesToMake;
	private double a = 8.82;
	private int numUnitCellsPerAxis = 5;
	private String initLatticeObjFileExtension = ".init.obj";
	
	/* **************************** */
	/* MAKE INITIAL LATTICE METHODS */
	/* **************************** */
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

		Lattice l = new Lattice(numUnitCellsPerAxis, a);
		l.setRandom(r);
		IdealTetrahedron[][][] lattice = null;
		InitialLattice init;
		MyPrintStream mps;
		writeToLog("\tcalling l.makeTetragonalFCCLattice() " + numLatticesToMake + " times.");
		
		for(int i = 0; i < numLatticesToMake; i++) {
			writeToLog("\tLattice " + (i+1) + " of " + numLatticesToMake);
			
			int numCycles = 1;
			
			switch(initLattice) {
			case RANDOM_SIXFOLD:
				lattice = l.makeTetragonalFCCLattice();
				writeToLog("\t\t made tetragonal (six options) FCC Lattice");
				break;
			case RANDOM_TOTALLY:
				lattice = l.makeRandomlyOrientedFCCLattice();
				writeToLog("\t\t made randomly oriented FCC Lattice");
				break;
			case MONOCLINIC_FIRST_SHELL:
				try {
					l.readBoxes(new File("sorted_eight.lattice"), Lattice.FIRST_SHELLS);
				} catch (FileNotFoundException e1) {
					e1.printStackTrace();
				}
				lattice = l.makeFirstShellLatticeRandomly(numCycles);
				writeToLog("\t\t made monoclinic first shell FCC Lattice with " + numCycles);
				break;
			case MONOCLINIC_SECOND_SHELL:
				try {
					File f = new File(".");
					System.out.println(f.getAbsolutePath());
					l.readBoxes(new File("sorted_rotated_second_shells.lattice"), Lattice.SECOND_SHELLS);
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
				lattice = l.makeSecondShellLatticeRandomly(1);
				writeToLog("\t\t made monoclinic second shell FCC Lattice with " + numCycles);
				break;
			case MONOCLINIC_STAR_BUILD:
				try {
					l.readBoxes(new File("110_stars.lattice"), Lattice.STARS);
				} catch (FileNotFoundException e1) {
					e1.printStackTrace();
				}
				double latticeEnergy = 0;
				int infiniteLoopEscapeMax = 25;
				int loopIdx = 0;
				lattice = l.makeStarLatticeRandomly(1);
				double zeroPotentialDistance = Z1Z2_distance / Math.pow(2., (1./6.));
				LennardJonesPotential ljp = new LennardJonesPotential(Z1, Z2, zeroPotentialDistance, getWellDepth());
				PotentialLookup ljl = new PotentialLookup(getLookupPrecision(), ljp);
				Simulate simul = new Simulate(ljl, getMinimumIntermolecularDistance(), getA());
				simul.setMaxTorque(getMaximumRotationalAnglePerRotation());
				simul.setMaxTrans(getMaximumTranslationalMotionPerMove());
				simul.setLattice3(lattice);
				simul.setRandom(new Random(getRandomNumberSeed()));
				simul.calculateTotalEnergy();
				latticeEnergy = simul.getTotalEnergy();
				while(true) {
					if(!Double.isInfinite(latticeEnergy))
						break;
					l.setRandom(new Random());
					lattice = l.makeStarLatticeRandomly(1);
					simul.setLattice3(lattice);
					simul.calculateTotalEnergy();
					latticeEnergy = simul.getTotalEnergy();
				}
				writeToLog("\t\t made monoclinic star FCC Lattice with " + numCycles);
			}
	
			
			init = new InitialLattice();
			init.setA(a);
			init.setLattice(lattice);
			init.setNumUnitCellsPerAxis(numUnitCellsPerAxis);
			init.setInitLatticeObjFileExtension(initLatticeObjFileExtension);
			
			writeObject(init, new File(objOutputFolder + File.separator + i + initLatticeObjFileExtension));
			
			if(outputAllTo != null)
				mps = new MyPrintStream(new File(outputAllTo + File.separator + ais.toString() + "--" + i + "--init.xyz"));
			else
				mps = new MyPrintStream(new File(xyzOutputFolder + File.separator + i + ".init.xyz"));
			LatticeTools.printLatticeXYZ(lattice, mps.getPrintStream(), a);
			mps.close();
			writeToLog("\t\tLattice written to file.");
			mps_log.flush();
		}
		writeToLog("Making " + numLatticesToMake + " initial lattices complete.");
	}
	/* ******************************** */
	/* END MAKE INITIAL LATTICE METHODS */
	/* ******************************** */
	
	/* *************************** */
	/* LATTICE PARAMETERS */
	/* *************************** */
	private ModifyLatticeTypes modLattice;
	private int startModifyIdx = 0;
	private int finishModifyIdx;
	private String modLatticeObjFileExtension = ".mod.obj";
	private double minimumIntermolecularDistance = 0;
	
	/* ********************** */
	/* interaction parameters */
	/* ********************** */
	private PotentialLookup ljl;
	private double lookupPrecision = 0.001;
	private int Z1 = 35;
	private int Z2 = 35;
	private double Z1Z2_distance = 3.8;
	private double zeroPotentialDistance;
	private double wellDepth = 0.034891;
	private String energyUnits = "eV";
	/* ********************** */
	/* random walk parameters */
	/* ********************** */
	private int numRandomWalksPerMolecule = 5;
	private int numStepsPerWalk = 256;
	private double maximumTranslationalMotionPerMove = 0.25;
	private double maximumRotationalAnglePerRotation = .25;
	private double energyCutoff = 0;
	private boolean useEnergyCutoff = false;
	private boolean disorderBetweenCycles = false;
	/* ******************************************** */
	/* monte carlo orientational changes parameters */
	/* ******************************************** */
	private int numTestsPerMolecule = 1;
	private boolean maxNumChangesEnabled = false;
	private int maxNumOrientationalChanges = 100000;
	
	/* ************************ */
	/* MODIFIED LATTICE METHODS */
	/* ************************ */
	private void setModLatticeParams(ModifiedLattice mod, int threadIdx) {
		mod.setModLattice(modLattice);
		mod.setModifiedLatticeObjFileName(modLatticeObjFileExtension);
		mod.setMinimumIntermolecularDistance(minimumIntermolecularDistance);
		
		/* interaction parameters */
		mod.setLookupPrecision(lookupPrecision);
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
		
		mod.setRandomNumberSeed(System.currentTimeMillis());
		
		mod.setInputLatticeObjFileName(objFiles[threadIdx].getAbsolutePath());
		mod.setOutputLatticeObjFileName(objOutputFolder + File.separator + threadIdx + modLatticeObjFileExtension);
		writeToLog("Thread: " + threadIdx + " object input file: " + mod.getInputLatticeObjFileName());
		if(outputAllTo != null)
			mod.setOutputLatticeXYZFileName(outputAllTo.getAbsolutePath()+ File.separator + ais.toString() + "--" + threadIdx + "--mod.xyz");
		else
			mod.setOutputLatticeXYZFileName(xyzOutputFolder + File.separator + threadIdx + ".xyz");
		mod.setThreadIdx(threadIdx);
	}
	/* **************************** */
	/* END MODIFIED LATTICE METHODS */
	/* **************************** */
	
	/* ************************** */
	/* ANALYZE LATTICE PARAMETERS */
	/* ************************** */
	
	/* *********************** */
	/* ANALYZE LATTICE METHODS */
	/* *********************** */
	private void analyzeLattice() {
		
	}
	/* *************************** */
	/* END ANALYZE LATTICE METHODS */
	/* *************************** */
	
	/* ******************************** */
	/* CALCULATE DIFFRACTION PARAMETERS */
	/* ******************************** */
	enum Projections {
		_1d,
		_2d_custom,
		_2d_001,
		_2d_001_90,
		_2d_010,
		_2d_100,
		_2d_011,
		_2d_111,
		_2d_211,
		_2d_210,
		_3d,
		;
	}
	enum DiffractionPrintStyle {
		xyI_column,
		array,
		;
	}
	private Projections diffCalcType;
	private DiffractionPrintStyle diffractionPrintStyle;
	private JVector[][] axes;
	private double qMaxX;
	private double qMaxY;
	private double qStep;
	private double zShift;
	private JPixel[] pixels;
	private int[] elemTypes;
	private double wavelength;
	private Lattices diffractionLattice;
	private double[][] xyI;
	private Vector<Projections> projections;
	private boolean areAtomsInCartesianCoordinates;
	
	/* ***************************** */
	/* CALCULATE DIFFRACTION METHODS */
	/* ***************************** */
	
	private void initMasterPixels() throws ClassNotFoundException, IOException {
		if(qStep == 0) {
			Lattices firstLattice = LatticeTools.readInLattice(objFiles[0]);
			qStep = 1. / ((double) firstLattice.getNumUnitCellsPerAxis());
			firstLattice = null;
		}
		pixels = DiffractionTools.initPixels_old(axes, qMaxX, qMaxY, qStep, elemTypes, wavelength, zShift);
	}
	private void setCalcDiffractionParams(CalculateDiffraction calc) {
		calc.setDiffCalcType(diffCalcType);
		calc.setDiffractionPrintStyle(diffractionPrintStyle);
		calc.setAxes(axes);
		calc.setqMaxX(qMaxX);
		calc.setqMaxY(qMaxY);
		calc.setqStep(qStep);
		calc.setPixels(pixels);
		calc.setElemTypes(elemTypes);
		calc.setWavelength(wavelength);
		calc.setXyI(xyI);
		calc.setareAtomsInCartesianCoordinates(areAtomsInCartesianCoordinates);
		if(diffractionOutputFolder == null) {
			File diffractionOutputFolder = new File(inputFolder.getAbsoluteFile().getParentFile() + File.separator + "diffraction");
			diffractionOutputFolder.mkdirs();
		}
		calc.setDiffractionOutputFolder(diffractionOutputFolder);
	}
	/*private void calcDiffraction(int i) throws ClassNotFoundException, IOException {
		Date startCalc = new Date();
		writeToLog("\t\tAbout to read in modified lattice: " + objFiles[i].getAbsolutePath());
		diffractionLattice = LatticeTools.readInLattice(objFiles[i]);			
		writeToLog("\t\tRead in modified lattice: " + objFiles[i].getAbsolutePath());
		Date endCalc = new Date();
		System.out.println("\t\tTotal time to read lattice: " + 
				(endCalc.getTime() - startCalc.getTime()) / 1000 + " s.");
		JAtom[] lattice = LatticeTools.molToFractionalAtoms(diffractionLattice.getLattice(), diffractionLattice.getA());
		writeToLog("\t\tConverted molecules to atoms.");
		startCalc = new Date();
		JavaToC.calcDiffraction(lattice, pixels, elemTypes);
		endCalc = new Date();
		writeToLog("\t\tCalculated diffraction of " + lattice.length + " atoms onto " + pixels.length + 
				" pixels. Total time: " + (endCalc.getTime() - startCalc.getTime()) / 1000 + " s.");
	}*/
	private void writeDiffractionToFile(CalculateDiffraction calc) {
		boolean appendToExistingFile = false;
		String outputFileRoot = calc.getObjectFileName();
		outputFileRoot = outputFileRoot.substring(outputFileRoot.lastIndexOf(File.separator)+1);
		File root = calc.getDiffractionOutputFolder();
		if(calc.getIndex() != null)
			root = new File(root + File.separator + calc.getIndex());
		root.mkdirs();
		MyPrintStream mps = new MyPrintStream(new File(root + File.separator + outputFileRoot + calc.getDiffCalcType().name() + ".xray"), appendToExistingFile);
		String str_axes = ((int) Math.rint(calc.getAxes()[0][0].i)) + "" + 
				((int) Math.rint(calc.getAxes()[0][0].j)) + "" + ((int) Math.rint(calc.getAxes()[0][0].k));
		mps.printCommentLine("Diffraction data for " + str_axes + " projection family");
		mps.printCommentLine(15);
		Lattices lattice = calc.getLattice();
		if(lattice != null) {
			mps.printCommentLine(15);
			mps.printCommentLine("Lattice Parameters:");
			String[] params = lattice.getParamsStringArray();
			for(String str : params)
				mps.printCommentLine(str);
			mps.printCommentLine(15);
			mps.printCommentLine(15);
		} else
			mps.printCommentLine("Diffraction calculated from XYZ file.");
		mps.printCommentLine(15);
		mps.printCommentLine("Diffraction Calculation Parameters: ");
		String[] params = calc.getDiffractionParamsStringArray();
		for(String str : params)
			mps.printCommentLine(str);

		mps.printCommentLine(15);
		mps.printCommentLine(15);
		mps.printCommentLine("Diffraction Calculation data");
		mps.printCommentLine(15);
		mps.printCommentLine("Diffraction Print Style: ");
		mps.println(calc.getDiffractionPrintStyle().name());
		switch(calc.getDiffractionPrintStyle()) {
		case array:
//			mps.println("X_MIN " + -1 * qMaxX);
//			mps.println("X_MAX " + qMaxX);
//			mps.println("Y_MIN " + -1 * qMaxY);
//			mps.println("Y_MAX " + qMaxY);
//			mps.println("Q_STEP " + qStep);
			writeToLog("Diffraction print style " + calc.getDiffractionPrintStyle().name() + 
				" is not yet enabled. Defaulting to " + DiffractionPrintStyle.xyI_column.name());
		case xyI_column:
			calc.pixelsToDoubleArray();
			String lines = "";
			double[][] xyI = calc.getXyI();
			for(int j = 0; j < calc.getXyI().length; j++) {
				if(j%1000 == 0) {
					mps.print(lines);
					lines = "";
				}
				switch(xyI[0].length) {
				case 2: 
					lines += "\n" + xyI[j][0] + "\t" + xyI[j][1];
					break;
				case 3:
					lines += "\n" + xyI[j][0] + "\t" + xyI[j][1] + "\t" + xyI[j][2];
					break;
				}
			}
			mps.println(lines);
			break;
		}
		mps.close();
	}
	private void initMasterXYI() {
		int numPixels = 0;
		for(double qy = -qMaxY; qy <= qMaxY; qy += qStep) {
			for(double qx = -qMaxX; qx <= qMaxX; qx += qStep) {
				numPixels++;
			}
		}
		xyI = new double[numPixels][3];
		int pixIdx = 0;
		for(double qy = -qMaxY; qy <= qMaxY; qy += qStep) {
			for(double qx = -qMaxX; qx <= qMaxX; qx += qStep, pixIdx++) {
				xyI[pixIdx][0] = Math.rint(qx / qStep);
				xyI[pixIdx][1] = Math.rint(qy / qStep);
				xyI[pixIdx][2] = 0;
			}
		}
		System.out.println("xyI done");
	}
	private void resetPixels() {
		for(int i = 0; i < pixels.length; i++)
			pixels[i].setI(0);
		
		for(int i = 0; i < xyI.length; i++)
			xyI[i][xyI[i].length-1] = 0;
	}
	/* ********************************* */
	/* END CALCULATE DIFFRACTION METHODS */
	/* ********************************* */
	
	
	private void init() {
		makeOutputFolder(outputRoot);
		initStreams();
		initInputFileArrays();
		r = new Random(randomNumberSeed);
		writeToLog("Random number seed: " + randomNumberSeed);
		writeToLog("Using: " + numThreadsToUse + " threads.");
	}
	
	private void finish() {
		closeStreams();
		pixels = null;
		xyI = null;
	}
	public Vector<Future<?>> run() {
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ssss");
		Date start = new Date();
		Vector<Future<?>> futures = new Vector<Future<?>>();
		switch(simulType) {
		case INITIAL_LATTICE:
			init();
			initObjOut();
			initXYZOut();
			if(numLatticesToMake == 0)
				numLatticesToMake = finishModifyIdx - startModifyIdx;
			writeToLog("RunSimulation.run() started at: " + dateFormat.format(start));
			makeInitialLattice();
			break;
		case MODIFY_LATTICE:
			init();
			initObjOut();
			initXYZOut();
			writeToLog("RunSimulation.run() started at: " + dateFormat.format(start));
			writeToLog("Starting thread idx: " + startModifyIdx);
			for(int i = startModifyIdx; i < finishModifyIdx; i++) {
				ModifiedLattice mod = new ModifiedLattice();
				setModLatticeParams(mod, i);
				Thread_ModifyLattice tml = new Thread_ModifyLattice(mod, mps_log);
				tml.addObserver(this);
				futures.add(es_simul.submit(tml));
				//es.execute(tml);
				writeToLog("Added thread " + i + " of " + finishModifyIdx);
			}
			break;
		case CALC_DIFFRACTION_OBJECT:
			init();
			writeToLog("RunSimulation.run() started at: " + dateFormat.format(start));
			try {
				initMasterPixels();
			} catch (ClassNotFoundException e1) {
				String errorMsg = "Could not load object file in to memory.";
				errorMsg += e1.getMessage();
				abruptExit(errorMsg);
			} catch (IOException e1) {
				String errorMsg = "Could not load object file in to memory.";
				errorMsg += e1.getMessage();
				abruptExit(errorMsg);
			}
			initMasterXYI();
			writeToLog("Starting image calc idx: " + startModifyIdx);
			for(int i = startModifyIdx; i < finishModifyIdx; i++) {
				String objectFileName = objFiles[i].getAbsolutePath() + "";
				CalculateDiffraction calc = new CalculateDiffraction(objectFileName);
				setCalcDiffractionParams(calc);
				calc.setIndex(ais.getName());
				Thread_CalculateDiffraction_Object calcDiff = new Thread_CalculateDiffraction_Object(calc, mps_log);
				calcDiff.addObserver(this);
				futures.add(es_diffraction.submit(calcDiff));
			}
			
			break;
		case ANALYZE_LATTICE:
			init();
			writeToLog("RunSimulation.run() started at: " + dateFormat.format(start));
			analyzeLattice();
			break;
		case CALC_DIFFRACTION_XYZ:
			init();
			initStreams();
			InitialLattice lattice = new InitialLattice();
			lattice.setA(a);
			lattice.setNumUnitCellsPerAxis(numUnitCellsPerAxis);
			writeToLog("RunSimulation.run() started at: " + dateFormat.format(start));
			try {
				initMasterPixels();
			} catch (ClassNotFoundException e1) {
				String errorMsg = "Could not load object file in to memory.";
				errorMsg += e1.getMessage();
				abruptExit(errorMsg);
			} catch (IOException e1) {
				String errorMsg = "Could not load object file in to memory.";
				errorMsg += e1.getMessage();
				abruptExit(errorMsg);
			}
			initMasterXYI();
			writeToLog("Starting image calc idx: " + startModifyIdx);
			for(int i = startModifyIdx; i < finishModifyIdx; i++) {
				if(!xyzFiles[i].isDirectory()) {
					String objectFileName = xyzFiles[i].getAbsolutePath() + "";
					CalculateDiffraction calc = new CalculateDiffraction(objectFileName);
					calc.setLattice(lattice);
					setCalcDiffractionParams(calc);
					Thread_CalculateDiffraction_XYZ calcDiff = new Thread_CalculateDiffraction_XYZ(calc, mps_log);
					calcDiff.addObserver(this);
					futures.add(es_diffraction.submit(calcDiff));
				}
			}
			break;
		default:
			break;
		}
		Date finish = new Date();
		writeToLog("RunSimulation.run() complete at: " + dateFormat.format(new Date()));
		writeToLog("Total time: " + Math.rint((finish.getTime() - start.getTime()) / 60000.*1000)/1000 + " minutes");
		return futures;
	}
	@Override
	public void update(Observable arg0, Object arg1) {
		if(arg1 instanceof String) {
			String str = (String) arg1;
			if(arg0 instanceof Thread_CalculateDiffraction_Object) {
				if(str.compareTo(Thread_CalculateDiffraction_Object.COMPLETE) == 0) {
					Thread_CalculateDiffraction_Object thread_calc = (Thread_CalculateDiffraction_Object) arg0;
					CalculateDiffraction calc = thread_calc.getCalc();
					writeDiffractionToFile(calc);
					thread_calc.getMpsLog().println(calc.getShortObjectFileName() + ": Diffraction written to file.");
					resetPixels();
					thread_calc.getMpsLog().println(calc.getShortObjectFileName() + ": Pixels reset.");
				} else if(str.compareTo(Thread_ModifyLattice.COMPLETE) == 0) {
					
				}
			} else if(arg0 instanceof Thread_CalculateDiffraction_XYZ) {
				if(str.compareTo(Thread_CalculateDiffraction_XYZ.COMPLETE) == 0) {
					Thread_CalculateDiffraction_XYZ thread_calc = (Thread_CalculateDiffraction_XYZ) arg0;
					CalculateDiffraction calc = thread_calc.getCalc();
					writeDiffractionToFile(calc);
					thread_calc.getMpsLog().println(calc.getShortObjectFileName() + ": Diffraction written to file.");
					resetPixels();
					thread_calc.getMpsLog().println(calc.getShortObjectFileName() + ": Pixels reset.");
					calc.nullify();
					thread_calc.nullify();
					thread_calc = null;
				}
			}
		} else if(arg0 instanceof WatchForNewFile) { 
			if(arg1 instanceof Object[]) {
				Object[] obj = (Object[]) arg1;
				String str = (String) obj[0];
				if(str.compareTo(WatchForNewFile.NEW_FILES_TO_COMPUTE) == 0) {
					if(obj.length > 1) {
						File[] files = (File[]) obj[1];
						inputFolder = files[0].getParentFile();
						xyzFiles = files;
						for(int j = 0; j < xyzFiles.length; j++) {
							startModifyIdx = j;
							finishModifyIdx = (j+1);
							runDiffractionCalc(this, SimulationToRun.CALC_DIFFRACTION_XYZ);
						}
					}
				}
			}
		} else {
			writeToLog(arg1.toString());
		}
	}
	
	private void abruptExit(String terminationMessage) {
		writeToLog("Simulation exiting: " + terminationMessage);
		System.exit(1);
	}
	private static File runRandomWalk(RunSimulationThread simul) {
		simul.simulType = SimulationToRun.MODIFY_LATTICE;

		simul.modLattice = ModifyLatticeTypes.RANDOM_WALK;
		simul.numRandomWalksPerMolecule = 1;
		//simul.run();

		block(simul.run());
		//simul.run();
		
		return simul.outputFolder;
	}
	
	private static File runMonteCarloOrientational(RunSimulationThread simul) {
		simul.simulType = SimulationToRun.MODIFY_LATTICE;

		simul.modLattice = ModifyLatticeTypes.MONTE_CARLO_ORIENTATIONAL;
		
		simul.numTestsPerMolecule = simul.numTestsPerMolecule;
		simul.maxNumChangesEnabled = simul.maxNumChangesEnabled;
		simul.maxNumOrientationalChanges = simul.maxNumOrientationalChanges;

		block(simul.run());
		
		return simul.outputFolder;
	}
	
	private static File makeInitialLattices(RunSimulationThread simul) {
		
		simul.simulType = SimulationToRun.INITIAL_LATTICE;
		block(simul.run());
		
		return simul.outputFolder;
	}
	private static void runDiffractionCalc(RunSimulationThread simul, SimulationToRun simulType) {
		for(Projections proj : simul.projections) {
			simul.simulType = simulType;
			simul.diffCalcType = proj;
			JVector vx = JVector.x;
			JVector vy = JVector.y;
			switch(simul.diffCalcType) {
			case _1d:
				System.out.println(simul.diffCalcType.name() + " is not yet set up.  Exiting calculation...");
				break;
			case _2d_001:
				vx = JVector.x;
				vy = JVector.y;
				simul.axes = JVector.axes100U_ZXY;
				break;
			case _2d_011:
				vx = new JVector(-1, 0, 0);
				vy = new JVector(0, 1, -1);
				simul.axes = JVector.axes110U_ZXY;
//				simul.axes = new JVector[][] {{vx, vy, JVector.cross(vx, vy)}};
//				simul.axes = new JVector[][] {{JVector.x, JVector.y, JVector.z}};
				break;
			case _2d_111:
				vx = new JVector(1, -1, 0);
				vy = new JVector(-1, -1, 2);
				simul.axes = JVector.axes111U_ZXY;
//				simul.axes = new JVector[][] {JVector.axes111_uniqueCrystallographic[0]};
//				simul.axes = new JVector[][] {{vx, vy, JVector.cross(vx, vy)}};
				break;
			case _2d_211:
				vx = new JVector(1, 1, 1);
				vy = new JVector(1, -1, 0);
//				simul.axes = JVector.getUAxes_XYZ(JVector.get111Family(), JVector.get110Family(), JVector.get112Family());;
				simul.axes = new JVector[][] {{vx, vy, JVector.cross(vx, vy)}};
				break;
			case _2d_210:
				vx = new JVector(0, 0, 1);
				vy = new JVector(2, 1, 0);
//				simul.axes = JVector.getUAxes_XYZ(JVector.get100Family(), JVector.get210Family(), JVector.get210Family());
				simul.axes = new JVector[][] {{vx, vy, JVector.cross(vx, vy)}};
				break;
			case _2d_custom:
				System.out.println(simul.diffCalcType.name() + " is not yet set up.  Exiting calculation...");
				System.exit(0);
			case _3d:
				System.out.println(simul.diffCalcType.name() + " is not yet set up.  Exiting calculation...");
				break;
			case _2d_010:
				vx = JVector.z;
				vy = JVector.x;
				simul.axes = new JVector[][] {{vx, vy, JVector.cross(vx, vy)}};
				break;
			case _2d_100:
				vx = JVector.y;
				vy = JVector.z;
				simul.axes = new JVector[][] {{vx, vy, JVector.cross(vx, vy)}};
				break;
			case _2d_001_90:
				vx = JVector.multiply(JVector.y, -1);
				vy = JVector.x;
				simul.axes = new JVector[][] {{vx, vy, JVector.cross(vx, vy)}};
				break;
			default:
				break;
			}

			simul.qMaxY = simul.qMaxX * vx.length() / vy.length();
			
			//Vector<Future<?>> futures = simul.run();
	//		System.out.println("About to call simul.run() for: " + simul.diffCalcType.name());
			block(simul.run());
			simul.finish();
	//		System.out.println("Block finished for: " + simul.diffCalcType.name());
//			simul.run();
		}
	}
	private static void block(Vector<Future<?>> futures) {
		int numRunning = futures.size();
		int millisToSleep = 1000;
		while(numRunning > 0) {
			numRunning = 0;
			for(Future<?> f : futures) {
				try {
					if( f.get() != null ) {
						break;
					}
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}
			}
			try {
				Thread.sleep(millisToSleep);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	public static void main(String[] args) {
		String drive = "";
		drive = "Z:\\";
		String root = drive + "Simulation" + File.separator + "Eric" + File.separator + "CBr4";
//		String root = "D:\\Documents referenced in lab notebooks\\";
		String user = "Dill";
		int pageNumber = 149;
		int notebookNumber = 4;
		String index = "b";
		
		File outputFolder = new File(root + File.separator + user + "-" + notebookNumber + File.separator + pageNumber);
		outputFolder.mkdirs();
		File inputFolder = new File(root + File.separator + "EDD_" + notebookNumber + "-" + pageNumber + index);
		long minutesToSleep = 0;
		try {
			Thread.sleep(minutesToSleep*60*1000);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		boolean monteCarloFirst = false;
		int walkLength = 256;
		int[] walkLengths = new int[] {walkLength, walkLength, walkLength, walkLength, 
				walkLength, walkLength, walkLength, walkLength, };//512, 512, 512, 512, 512, 512, 512, 512, 512, 512}; 
//		int[] walkLengths = new int[] {};
		inputFolder = null;
		RunSimulationThread simul = new RunSimulationThread(outputFolder, inputFolder);

		int startModifyIdx = 0;
		int finishModifyIdx = 100;
		simul.initLattice = InitialLatticeTypes.MONOCLINIC_STAR_BUILD;
		simul.numLatticesToMake = finishModifyIdx - startModifyIdx;
		simul.pageNumber = pageNumber;
		simul.notebookNumber = notebookNumber;
		simul.numThreadsToUse = 3;
		simul.numUnitCellsPerAxis = 10;
		simul.maximumTranslationalMotionPerMove = .1;
		
		simul.elemTypes = new int[] {6, 35};
		simul.diffractionPrintStyle = DiffractionPrintStyle.xyI_column;
		simul.wavelength = 0.13702;
		simul.initExecutors();
		Vector<File> inputFolders = new Vector<File>();
//		inputFolders.add(inputFolder);
		Vector<Projections> projections = new Vector<Projections>();
		projections.add(Projections._2d_001);
		projections.add(Projections._2d_011);
		projections.add(Projections._2d_111);
//		projections.add(Projections._2d_211);
//		projections.add(Projections._2d_210);
		simul.outputAllTo = new File(root + File.separator + user + "-" + notebookNumber + File.separator + pageNumber + File.separator + "simulOut");
		simul.outputAllTo.mkdirs();
//		simul.outputAllTo = new File("Z:\\Simulation\\Eric\\CBr4\\Diffraction patterns will be automatically calculated in this folder");

		/** UNCOMMENT THIS BLOCK TO CALCULATE DIFFRACTION PATTERNS IN A SPECIFIC SET OF FOLDERS */
		projections.add(Projections._2d_001);
		
		simul.projections = projections;
		simul.qMaxX = 10;
		simul.qStep = 1. / 10.;
		simul.a = 	1;
		simul.zShift = 2;
		double maxShift = 2.7;
		double shiftStep = 2*simul.qStep;
		double startShift = 2.6;
		
		int numToCalculate = (int) Math.rint((maxShift - startShift) / shiftStep);
		
		inputFolder = new File("D:\\Documents referenced in lab notebooks\\Dill-4\\149\\b");
		inputFolders.add(inputFolder);
		for(int idx = 0; idx < numToCalculate; idx++) {
			double shift = idx * shiftStep + startShift;
			if(shift > 0 && shift < 1.4)
				continue;
			shift = Math.rint(shift*100)/100.;
			simul.zShift = shift;
			
			int step = 1;
			simul.diffractionOutputFolder = new File(inputFolder.getParentFile() + File.separator + "diffraction -- shifted up " + shift);
			simul.outputFolder = simul.diffractionOutputFolder;
			
			simul.diffractionOutputFolder.mkdirs();
			simul.outputFolder.mkdirs();
			startModifyIdx = 0;
			for(File inFolder : inputFolders) {
				simul.inputFolder = inFolder;
				finishModifyIdx = new File(inFolder + File.separator + "xyz").listFiles().length;
				for(int j = startModifyIdx; j < finishModifyIdx; j+= step) {
					simul.startModifyIdx = j;
					simul.finishModifyIdx = j+step;
					runDiffractionCalc(simul, SimulationToRun.CALC_DIFFRACTION_XYZ);
				}
			}
		}

		/** UNCOMMENT THIS BLOCK TO CALCULATE DIFFRACTION PATTERNS AS OUTPUT FILES ARE DUMPED INTO A FOLDER */
		
//		File folderToWatch = simul.outputAllTo;
//		folderToWatch = new File("D:\\Documents referenced in lab notebooks\\Dill-4\\146\\simulOut");
//		simul.areAtomsInCartesianCoordinates = false;
//		simul.a = 1;
//		simul.qMaxX = 10;
//		simul.projections = projections;
//		simul.qStep = 1. / 30.;
//		simul.outputRoot = folderToWatch.getParentFile();
//		simul.outputFolder = new File(simul.outputRoot + File.separator + "diffraction");
//		WatchForNewFile wf = new WatchForNewFile(folderToWatch);
//		wf.addObserver( simul);
		
		/** UNCOMMENT THIS BLOCK TO RUN SIMULATIONS WHICH WILL DUMP .XYZ FILES INTO A FOLDER */
		
//		File input = new File(drive + File.separator + "Simulation" + File.separator + "Eric" + File.separator + "CBr4" + File.separator + "Dill-4" + File.separator + "95" + File.separator + "EDD_4-95b");
//		inputFolders.add(input);
		
//		int monteCarloOffset = 0;
//		for(int i = 0; i < walkLengths.length+1; i++) {
//			simul.startModifyIdx = startModifyIdx;
//			simul.finishModifyIdx = finishModifyIdx;
//			if(i == 0) {
//				inputFolders.add(makeInitialLattices(simul));
//			} else {
//				if(monteCarloFirst) {
//					simul.inputFolder = inputFolders.get(i-1);
//					simul.numStepsPerWalk = walkLengths[i-1];
//					inputFolders.add(runMonteCarloOrientational(simul));
//					monteCarloFirst = false;
//					monteCarloOffset = 1;
//					i--;
//				} else {
//					simul.inputFolder = inputFolders.get(i-1 + monteCarloOffset);
//					simul.numStepsPerWalk = walkLengths[i-1];
//					inputFolders.add(runRandomWalk(simul));
//				}
//			}
//		}
		
		/** **************************** */
		
//		for(int i = 0; i < walkLengths.length+1; i++) {
//			simul.startModifyIdx = startModifyIdx;
//			simul.finishModifyIdx = finishModifyIdx;
//			if(i == 0) {
//				simul.initLattice = InitialLatticeTypes.MONOCLINIC_SECOND_SHELL;
////				simul.initLattice = InitialLatticeTypes.RANDOM_TOTALLY;
//				inputFolders.add(makeInitialLattices(simul));
//			} else {
//				simul.inputFolder = inputFolders.get(i-1);
//				simul.numStepsPerWalk = walkLengths[i-1];
//				inputFolders.add(runRandomWalk(simul));
//			}
//		}
//			int step = 5;
//			simul.inputFolder = inputFolders.get(i);
//			for(int j = startModifyIdx; j < finishModifyIdx; j+= step) {
//				simul.startModifyIdx = j;
//				simul.finishModifyIdx = j+step;
//				runDiffractionCalc(simul, SimulationToRun.CALC_DIFFRACTION_OBJECT);
//			}
//		}
//		
//		simul.ais = new AlphabeticIndexingSystem(index);
//		runRandomWalk(simul);
//		
////		int step = 5;
//		simul.inputFolder = new File("Simulation\\Eric\\CBr4\\Dill-4\\89\\obj\\xyzOut\\all other times");
//		simul.objFiles = simul.inputFolder.listFiles();
////		simul.outputFolder = new File(simul.inputFolder + File.separator + "xrayCalc");
//		simul.numUnitCellsPerAxis = 15;
//		simul.qStep = 1. / (double) simul.numUnitCellsPerAxis;
//		startModifyIdx = 0;
//		finishModifyIdx = simul.objFiles.length;
//		int numSteps = (int) (finishModifyIdx / step);
//		int j = 0;
//		for(j = 0; j < numSteps; j++) {
//			simul.startModifyIdx = j*step;
//			simul.finishModifyIdx = (j+1)*+step;
//			runDiffractionCalc(simul, SimulationToRun.CALC_DIFFRACTION_XYZ);
//		}
//		simul.startModifyIdx = j;
//		simul.finishModifyIdx = finishModifyIdx;
//		runDiffractionCalc(simul, SimulationToRun.CALC_DIFFRACTION_XYZ);

	}
	
	/* ******************* */
	/* GETTERS AND SETTERS */
	/* ******************* */
	public File[] getObjFiles() { return objFiles; }
	public void setObjFiles(File[] objFiles) { this.objFiles = objFiles; }
	public File getInputFolder() { return inputFolder; }
	public void setInputFolder(File inputFolder) { this.inputFolder = inputFolder; }
	public File getOutputRoot() { return outputRoot; }
	public void setOutputRoot(File outputRoot) { this.outputRoot = outputRoot; }
	public int getNumThreadsToUse() { return numThreadsToUse; }
	public void setNumThreadsToUse(int numThreadsToUse) { this.numThreadsToUse = numThreadsToUse; }
	public SimulationToRun getSimulType() { return simulType; }
	public void setSimulType(SimulationToRun simulType) { this.simulType = simulType; }
	public Random getR() { return r; }
	public void setR(Random r) { this.r = r; }
	public long getRandomNumberSeed() { return randomNumberSeed; }
	public void setRandomNumberSeed(long randomNumberSeed) { this.randomNumberSeed = randomNumberSeed; }
	public ExecutorService getEs_simul() { return es_simul; }
	public void setEs_simul(ExecutorService es_simul) { this.es_simul = es_simul; }
	public ExecutorService getEs_diffraction() { return es_diffraction; }
	public void setEs_diffraction(ExecutorService es_diffraction) { this.es_diffraction = es_diffraction;}
	public DateFormat getDf() {  return df; }
	public void setDf(DateFormat df) { this.df = df; }
	public File getOutputFolder() { return outputFolder; }
	public void setOutputFolder(File outputFolder) { this.outputFolder = outputFolder; }
	public File getObjOutputFolder() { return objOutputFolder; }
	public void setObjOutputFolder(File objOutputFolder) { this.objOutputFolder = objOutputFolder; }
	public File getAnalysisOutputFolder() { return analysisOutputFolder; }
	public void setAnalysisOutputFolder(File analysisOutputFolder) { this.analysisOutputFolder = analysisOutputFolder; }
	public File getDiffractionOutputFolder() { return diffractionOutputFolder; }
	public void setDiffractionOutputFolder(File diffractionOutputFolder) { this.diffractionOutputFolder = diffractionOutputFolder; }
	public MyPrintStream getMps_log() { return mps_log; }
	public void setMps_log(MyPrintStream mps_log) { this.mps_log = mps_log; }
	public MyObjectOutputStream getMoos() { return moos; }
	public void setMoos(MyObjectOutputStream moos) { this.moos = moos; }
	public String getInitials() { return initials; }
	public void setInitials(String initials) { this.initials = initials; }
	public int getNotebookNumber() { return notebookNumber; }
	public void setNotebookNumber(int notebookNumber) { this.notebookNumber = notebookNumber; }
	public int getPageNumber() { return pageNumber; }
	public void setPageNumber(int pageNumber) { this.pageNumber = pageNumber; }
	public AlphabeticIndexingSystem getAis() { return ais; }
	public void setAis(AlphabeticIndexingSystem ais) { this.ais = ais; }
	public InitialLatticeTypes getInitLattice() { return initLattice; }
	public void setInitLattice(InitialLatticeTypes initLattice) { this.initLattice = initLattice; }
	public int getNumLatticesToMake() { return numLatticesToMake; }
	public void setNumLatticesToMake(int numLatticesToMake) { this.numLatticesToMake = numLatticesToMake;}
	public double getA() { return a; }
	public void setA(double a) { this.a = a; }
	public int getNumUnitCellsPerAxis() { return numUnitCellsPerAxis; }
	public void setNumUnitCellsPerAxis(int numUnitCellsPerAxis) { this.numUnitCellsPerAxis = numUnitCellsPerAxis; }
	public String getInitLatticeObjFileExtension() { return initLatticeObjFileExtension; }
	public void setInitLatticeObjFileExtension(String initLatticeObjFileExtension) { this.initLatticeObjFileExtension = initLatticeObjFileExtension; }
	public ModifyLatticeTypes getModLattice() { return modLattice; }
	public void setModLattice(ModifyLatticeTypes modLattice) { this.modLattice = modLattice; }
	public int getStartModifyIdx() { return startModifyIdx; }
	public void setStartModifyIdx(int startModifyIdx) { this.startModifyIdx = startModifyIdx; }
	public int getFinishModifyIdx() { return finishModifyIdx; }
	public void setFinishModifyIdx(int finishModifyIdx) { this.finishModifyIdx = finishModifyIdx; }
	public String getModLatticeObjFileExtension() { return modLatticeObjFileExtension; }
	public void setModLatticeObjFileExtension(String modLatticeObjFileExtension) { this.modLatticeObjFileExtension = modLatticeObjFileExtension; }
	public double getMinimumIntermolecularDistance() { return minimumIntermolecularDistance; }
	public void setMinimumIntermolecularDistance( double minimumIntermolecularDistance) { this.minimumIntermolecularDistance = minimumIntermolecularDistance; }
	public double getLookupPrecision() { return lookupPrecision; }
	public void setLookupPrecision(double lookupPrecision) { this.lookupPrecision = lookupPrecision; }
	public double getWellDepth() { return wellDepth; }
	public void setWellDepth(double wellDepth) { this.wellDepth = wellDepth; }
	public String getEnergyUnits() { return energyUnits; }
	public void setEnergyUnits(String energyUnits) { this.energyUnits = energyUnits; }
	public int getNumRandomWalksPerMolecule() { return numRandomWalksPerMolecule; }
	public void setNumRandomWalksPerMolecule(int numRandomWalksPerMolecule) { this.numRandomWalksPerMolecule = numRandomWalksPerMolecule; }
	public int getNumStepsPerWalk() { return numStepsPerWalk; }
	public void setNumStepsPerWalk(int numStepsPerWalk) { this.numStepsPerWalk = numStepsPerWalk; }
	public double getMaximumTranslationalMotionPerMove() { return maximumTranslationalMotionPerMove; }
	public void setMaximumTranslationalMotionPerMove(double maximumTranslationalMotionPerMove) { this.maximumTranslationalMotionPerMove = maximumTranslationalMotionPerMove; }
	public double getMaximumRotationalAnglePerRotation() { return maximumRotationalAnglePerRotation; }
	public void setMaximumRotationalAnglePerRotation(double maximumRotationalAnglePerRotation) { this.maximumRotationalAnglePerRotation = maximumRotationalAnglePerRotation; }
	public double getEnergyCutoff() { return energyCutoff; }
	public void setEnergyCutoff(double energyCutoff) { this.energyCutoff = energyCutoff; }
	public boolean isUseEnergyCutoff() { return useEnergyCutoff; }
	public void setUseEnergyCutoff(boolean useEnergyCutoff) { this.useEnergyCutoff = useEnergyCutoff; }
	public boolean isDisorderBetweenCycles() { return disorderBetweenCycles; }
	public void setDisorderBetweenCycles(boolean disorderBetweenCycles) { this.disorderBetweenCycles = disorderBetweenCycles; }
	public int getNumTestsPerMolecule() { return numTestsPerMolecule; }
	public void setNumTestsPerMolecule(int numTestsPerMolecule) { this.numTestsPerMolecule = numTestsPerMolecule; }
	public boolean isMaxNumChangesEnabled() { return maxNumChangesEnabled; }
	public void setMaxNumChangesEnabled(boolean maxNumChangesEnabled) { this.maxNumChangesEnabled = maxNumChangesEnabled; }
	public int getMaxNumOrientationalChanges() { return maxNumOrientationalChanges; }
	public void setMaxNumOrientationalChanges(int maxNumOrientationalChanges) { this.maxNumOrientationalChanges = maxNumOrientationalChanges; }
	public Projections getDiffCalcType() { return diffCalcType; }
	public void setDiffCalcType(Projections diffCalcType) { this.diffCalcType = diffCalcType; }
	public DiffractionPrintStyle getDiffractionPrintStyle() { return diffractionPrintStyle; }
	public void setDiffractionPrintStyle(DiffractionPrintStyle diffractionPrintStyle) { this.diffractionPrintStyle = diffractionPrintStyle; }
	public JVector[][] getAxes() { return axes; }
	public void setAxes(JVector[][] axes) { this.axes = axes; }
	public double getqMaxX() { return qMaxX; }
	public void setqMaxX(double qMaxX) { this.qMaxX = qMaxX; }
	public double getqMaxY() { return qMaxY; }
	public void setqMaxY(double qMaxY) { this.qMaxY = qMaxY; }
	public double getqStep() { return qStep; }
	public void setqStep(double qStep) { this.qStep = qStep; }
	public JPixel[] getPixels() { return pixels; }
	public void setPixels(JPixel[] pixels) { this.pixels = pixels; }
	public int[] getElemTypes() { return elemTypes; }
	public void setElemTypes(int[] elemTypes) { this.elemTypes = elemTypes; }
	public double getWavelength() { return wavelength; }
	public void setWavelength(double wavelength) { this.wavelength = wavelength; }
	public Lattices getDiffractionLattice() { return diffractionLattice; }
	public void setDiffractionLattice(Lattices diffractionLattice) { this.diffractionLattice = diffractionLattice; }
}