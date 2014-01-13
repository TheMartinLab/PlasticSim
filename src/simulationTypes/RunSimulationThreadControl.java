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
import io.StringConverter;

import java.io.File;
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
import chemistry.JAtomTools;
import defaultPackage.IdealTetrahedron;
import defaultPackage.JAtom;
import defaultPackage.JPixel;
import defaultPackage.JVector;
import defaultPackage.JavaToC;
import defaultPackage.Lattice;
import defaultPackage.PotentialLookup;

public class RunSimulationThreadControl implements Observer {

	/* ************ */
	/* CONSTRUCTORS */
	/* ************ */
	
	public RunSimulationThreadControl(File outputRoot) {
		this.outputRoot = outputRoot;
	}
	public RunSimulationThreadControl(File outputRoot, File inputFolder) {
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
	private SimulationToRun simulType;
	private Random r;
	private long randomNumberSeed = System.currentTimeMillis();
	private ExecutorService es_simul, es_diffraction;
	private DateFormat df = DateFormat.getInstance();
	enum SimulationToRun {
		INITIAL_LATTICE,
		MODIFY_LATTICE,
		CALC_DIFFRACTION,
		ANALYZE_LATTICE,
		;
	}
	
	/* ******************************** */
	/* INITIALIZATION AND INPUT METHODS */
	/* ******************************** */
	private void initExecutors() {
		es_simul = Executors.newFixedThreadPool(numThreadsToUse);
		es_diffraction = Executors.newFixedThreadPool(1);
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
	private void initStreams() {
		if(outputFolder == null)
			makeOutputFolder(outputRoot);
			
		outputFolder.mkdirs();
		
		mps_log = new MyPrintStream(new File(outputFolder + File.separator + getExptName() + ".log"));
		
		writeToLog("Working Directory: " + mps_log.getFile().getAbsolutePath());
		writeToLog("Simulation type: " + simulType.name());
	}
	private void makeOutputFolder(File outputRoot) {
		switch(simulType) {
		case ANALYZE_LATTICE:
		case CALC_DIFFRACTION:
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
		writeToLog("\t\tModified lattice object written to file.");
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
	private int numLatticesToMake = 1;
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
		switch(initLattice) {
		case RANDOM_SIXFOLD:
			Lattice l = new Lattice(numUnitCellsPerAxis, a);
			l.setRandom(r);
			IdealTetrahedron[][][] lattice;
			InitialLattice init;
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
				
				writeObject(init, new File(objOutputFolder + File.separator + i + initLatticeObjFileExtension));
				
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
	/* ******************************** */
	/* END MAKE INITIAL LATTICE METHODS */
	/* ******************************** */
	
	/* *************************** */
	/* MODIFIED LATTICE PARAMETERS */
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
	private double wellDepth = 0.01163;
	private String energyUnits = "eV";
	/* ********************** */
	/* random walk parameters */
	/* ********************** */
	private int numRandomWalksPerMolecule = 5;
	private int numStepsPerWalk = 256;
	private double maximumTranslationalMotionPerMove = 0.25;
	private double maximumRotationalAnglePerRotation = 1;
	private double energyCutoff = 0;
	private boolean useEnergyCutoff = false;
	private boolean disorderBetweenCycles = false;
	/* ******************************************** */
	/* monte carlo orientational changes parameters */
	/* ******************************************** */
	private int numTestsPerMolecule = 5;
	private boolean maxNumChangesEnabled = false;
	private int maxNumOrientationalChanges = 1000;
	
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
	public enum Projections {
		_1d,
		_2d_custom,
		_2d_001,
		_2d_011,
		_2d_111,
		_3d,
		;
	}
	public enum DiffractionPrintStyle {
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
	private JPixel[] pixels;
	private int[] elemTypes;
	private double wavelength;
	private Lattices diffractionLattice;
	private double[][] xyI;
	
	/* ***************************** */
	/* CALCULATE DIFFRACTION METHODS */
	/* ***************************** */
	
	private void initMasterPixels() throws ClassNotFoundException, IOException {
		Lattices firstLattice = LatticeTools.readInLattice(objFiles[0]);
		if(qStep == 0)
			qStep = 1. / ((double) firstLattice.getNumUnitCellsPerAxis());
		pixels = DiffractionTools.initPixels(axes, qMaxX, qMaxY, qStep, elemTypes, wavelength);
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
		File diffractionOutputFolder = new File(inputFolder.getAbsoluteFile() + File.separator + "diffraction");
		diffractionOutputFolder.mkdirs();
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
		
		MyPrintStream mps = new MyPrintStream(new File(calc.getDiffractionOutputFolder() + File.separator + outputFileRoot + calc.getDiffCalcType().name() + ".xray"), appendToExistingFile);
		String str_axes = ((int) Math.rint(calc.getAxes()[0][0].i)) + "" + 
				((int) Math.rint(calc.getAxes()[0][0].j)) + "" + ((int) Math.rint(calc.getAxes()[0][0].k));
		mps.printCommentLine("Diffraction data for " + str_axes + " projection family");
		mps.printCommentLine(15);
		mps.printCommentLine(15);
		mps.printCommentLine("Lattice Parameters:");
		String[] params = calc.getLattice().getParamsStringArray();
		for(String str : params)
			mps.printCommentLine(str);
		mps.printCommentLine(15);
		mps.printCommentLine(15);
		mps.printCommentLine(15);
		mps.printCommentLine("Diffraction Calculation Parameters: ");
		params = calc.getDiffractionParamsStringArray();
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
	
	
	public Vector<Future<?>> run() {
		makeOutputFolder(outputRoot);
		initStreams();
		initInputFileArrays();
		Date start = new Date();
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		writeToLog("RunSimulation.run() started at: " + dateFormat.format(start));
		r = new Random(randomNumberSeed);
		writeToLog("Random number seed: " + randomNumberSeed);
		writeToLog("Using: " + numThreadsToUse + " threads.");
		Vector<Future<?>> futures = new Vector<Future<?>>();
		switch(simulType) {
		case INITIAL_LATTICE:
			initObjOut();
			initXYZOut();
			numLatticesToMake = finishModifyIdx - startModifyIdx;
			makeInitialLattice();
			break;
		case MODIFY_LATTICE:
			initObjOut();
			initXYZOut();
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
		case CALC_DIFFRACTION:
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
				Thread_CalculateDiffraction calcDiff = new Thread_CalculateDiffraction(calc, mps_log);
				calcDiff.addObserver(this);
				futures.add(es_diffraction.submit(calcDiff));
			}
			
			break;
		case ANALYZE_LATTICE:
			analyzeLattice();
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
			if(str.compareTo(Thread_CalculateDiffraction.COMPLETE) == 0) {
				Thread_CalculateDiffraction thread_calc = (Thread_CalculateDiffraction) arg0;
				CalculateDiffraction calc = thread_calc.getCalc();
				writeDiffractionToFile(calc);
				thread_calc.getMpsLog().println(calc.getShortObjectFileName() + ": Diffraction written to file.");
				resetPixels();
				thread_calc.getMpsLog().println(calc.getShortObjectFileName() + ": Pixels reset.");
			} else if(str.compareTo(Thread_ModifyLattice.COMPLETE) == 0) {
				
			}
		} else {
			writeToLog(arg1.toString());
		}
	}
	
	private void abruptExit(String terminationMessage) {
		writeToLog("Simulation exiting: " + terminationMessage);
		System.exit(1);
	}
	private static File runRandomWalk(RunSimulationThreadControl simul) {
		simul.simulType = SimulationToRun.MODIFY_LATTICE;

		simul.modLattice = ModifyLatticeTypes.RANDOM_WALK;
		simul.numRandomWalksPerMolecule = 1;
		//simul.run();

		block(simul.run());
		//simul.run();
		
		return simul.outputFolder;
	}
	
	private static File makeInitialLattices(RunSimulationThreadControl simul) {
		
		simul.simulType = SimulationToRun.INITIAL_LATTICE;
		simul.initLattice = InitialLatticeTypes.RANDOM_SIXFOLD;
		block(simul.run());
		
		return simul.outputFolder;
	}
	private static void runDiffractionCalc(RunSimulationThreadControl simul, Vector<Projections> projections) {
		for(Projections proj : projections) {
			simul.qMaxX = 10;
			simul.qStep = 1. / (double) simul.numUnitCellsPerAxis;
			simul.simulType = SimulationToRun.CALC_DIFFRACTION;
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
				simul.axes = JVector.axes100U_XYZ;
				break;
			case _2d_011:
				vx = new JVector(-1, 0, 0);
				vy = new JVector(0, 1, -1);
				simul.axes = JVector.axes110U_XYZ;
				//simul.axes = new JVector[][] {{vx, vy, JVector.cross(vx, vy)}};
				break;
			case _2d_111:
				vx = new JVector(1, -1, 0);
				vy = new JVector(-1, -1, 2);
				simul.axes = JVector.axes111U_XYZ;
				//simul.axes = new JVector[][] {{vx, vy, JVector.cross(vx, vy)}};
				break;
			case _2d_custom:
				System.out.println(simul.diffCalcType.name() + " is not yet set up.  Exiting calculation...");
				System.exit(0);
			case _3d:
				System.out.println(simul.diffCalcType.name() + " is not yet set up.  Exiting calculation...");
				break;
			}

			simul.qMaxY = simul.qMaxX * vx.length() / vy.length();
			
			//Vector<Future<?>> futures = simul.run();
	//		System.out.println("About to call simul.run() for: " + simul.diffCalcType.name());
	//		block(simul.run());
	//		System.out.println("Block finished for: " + simul.diffCalcType.name());
			simul.run();
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
		String root = "Z:\\Simulation\\Eric\\CBr4\\";
		String user = "Dill";
		int pageNumber = 88;
		int notebookNumber = 4;
		String index = "a";
		
		File outputFolder = new File(root + user + "-" + notebookNumber + File.separator + pageNumber);
		long minutesToSleep = 0;
		try {
			Thread.sleep(minutesToSleep*60*1000);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		int[] walkLengths = new int[] {512, 512, 512, 512, 512}; 
		RunSimulationThreadControl simul = new RunSimulationThreadControl(outputFolder);

		int startModifyIdx = 0;
		int finishModifyIdx = 25;
		simul.pageNumber = pageNumber;
		simul.notebookNumber = notebookNumber;
		simul.numThreadsToUse = 2;
		simul.numUnitCellsPerAxis = 15;
		
		simul.elemTypes = new int[] {6, 35};
		simul.diffractionPrintStyle = DiffractionPrintStyle.xyI_column;
		simul.wavelength = 0.13702;
		simul.numUnitCellsPerAxis = simul.numUnitCellsPerAxis;
		simul.initExecutors();
		Vector<File> inputFolders = new Vector<File>();
		Vector<Projections> projections = new Vector<Projections>();
		projections.add(Projections._2d_001);
		projections.add(Projections._2d_011);
		projections.add(Projections._2d_111);

		for(int i = 0; i < walkLengths.length+1; i++) {
			simul.startModifyIdx = startModifyIdx;
			simul.finishModifyIdx = finishModifyIdx;
			if(i == 0) {
				inputFolders.add(makeInitialLattices(simul));
			} else {
				simul.inputFolder = inputFolders.get(i-1);
				simul.numStepsPerWalk = walkLengths[i-1];
				inputFolders.add(runRandomWalk(simul));
			}

			int step = 5;
			simul.inputFolder = inputFolders.get(i);
			for(int j = startModifyIdx; j < finishModifyIdx; j+= step) {
				simul.startModifyIdx = j;
				simul.finishModifyIdx = j+step;
				runDiffractionCalc(simul, projections);
			}
		}
		
		simul.es_simul.shutdown();
		simul.es_diffraction.shutdown();
		//simul.ais = new AlphabeticIndexingSystem(index);

		//runRandomWalk(simul);
	}

}