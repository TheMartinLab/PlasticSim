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
import java.io.IOException;
import java.util.Vector;
import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import defaultPackage.HookePotential;
import defaultPackage.LennardJonesPotential;
import defaultPackage.PotentialLookup;
import simulationTools.LatticeTools;

public class RunDistortTetrahedra {

	public void run() {
		
	}
	
	
	public static void main(String[] args) {
		
		/* ********************* */
		/* SET UP THE POTENTIALS */
		/* ********************* */
		/* ************************** */
		/* HOOKE POTENTIAL PARAMETERS */
		/* ************************** */
		double tolerance = 1e-4;
		
		HookePotential hookeCBrIntramolec, hookeBrBrIntramolec;
		double cBrDist = 1.91;
		double brBrIntramolecDist = 3.119;
		// nu = 267 cm-1 from handbook of chemistry and physics
		// mu = 12.01 * 79.9 / (12.01 + 79.9) / 1000 / 6.022e23 = 1.73e-26
		// k = (2*PI*3e10*267)^2 * 1.73e-26 = 43.82
		double cBrSpringConstant;
		cBrSpringConstant = 43.82;	// upper Limit
//		cBrSpringConstant /= 4;		// lower limit
		double brBrSpringConstant = cBrSpringConstant / 10.;
//		cBrSpringConstant /= 10.;
		
		/* ********************************** */
		/* LENNARD JONES POTENTIAL PARAMETERS */
		/* ********************************** */
		LennardJonesPotential ljp_intermolec;
		double lookupPrecision = 0.001;
		int Z1 = 35;
		int Z2 = 35;
		double Z1Z2_distance = 3.8;
		double zeroPotentialDistance;
		double wellDepth = 0.034891;
		wellDepth *= 10;
		
			/* INIT POTENTIALS */
		zeroPotentialDistance = Z1Z2_distance / Math.pow(2., (1./6.));
		ljp_intermolec = new LennardJonesPotential(Z1, Z2, zeroPotentialDistance, wellDepth);
		hookeCBrIntramolec = new HookePotential(6, 35, cBrDist, cBrSpringConstant);
		hookeBrBrIntramolec = new HookePotential(35, 35, brBrIntramolecDist, brBrSpringConstant);
		
		/* *************** */
		/* MISC PARAMETERS */
		/* *************** */
			/* TIME */
		double timeStep = 0.01;
		int stepsToCheckEnergy = 100;
			/* OUTPUT */
		boolean outputXYZs = true;
		boolean outputMovie = true;
		int outputXYZsEveryThisManyTimeSteps = 100;
		int outputMovieXYZsEveryThisManyTimeSteps = 100;
		File outputFolder = null;
		double intermolecularMaxDist = 10;
			/* END CONDITIONS */
		boolean isEndingBecauseOfSmallMovement = false;
		double smallEndMovement = .0001;
		boolean isEndingBecauseOfLargeMovement = true;
		double largeEndMovement = .1;
		boolean isEndingBecauseOfEnergy = false;
		double endDeltaE = .005;
		endDeltaE *= stepsToCheckEnergy;
		boolean isEndingBecauseOfTimeSteps = true;
		int endTimeStep = 5000;
		
		boolean inEclipse = args.length == 0;
		int numProcessorsToUse = 1;
		File inputFolder = new File("D:\\Documents referenced in lab notebooks\\Dill-4\\89\\EDD_4-89d\\obj");
		if(!inEclipse) { 
			inputFolder = new File("obj");
			System.out.println("Looking for object folder in the .jar directory: " + inputFolder.getAbsolutePath());
			if(inputFolder.isDirectory())
				System.out.println("Object folder found. Using this as input.");
			else if(args.length > 1) {
				inputFolder = new File(args[1]);
				System.out.println("Default folder not found. Using user input from command line: " + inputFolder.getAbsolutePath());
			} else if(!inputFolder.isDirectory()) {
				System.out.println("Command line input: " + inputFolder.getAbsolutePath() + " is not a valid directory.");
				System.out.println("\n\nExiting...");
				System.exit(0);
			}
			numProcessorsToUse = Integer.valueOf(args[0]);
		}

		File[] inputObjects = inputFolder.listFiles();
		ExecutorService es = Executors.newFixedThreadPool(numProcessorsToUse);
		System.out.println("Using " + numProcessorsToUse + " processors.");
		Vector<Future<?>> futures = new Vector<Future<?>>();
		for(File f : inputObjects) {
			if(f.isDirectory())
				continue;
			
			String name = f.getName();
			name.substring(0, name.indexOf("."));
			/* ******************* */
			/* READ IN THE LATTICE */
			/* ******************* */
			
			Lattices lattice = null;
			try {
				lattice = LatticeTools.readInLattice(f);
			} catch (ClassNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			/* *************** */
			/* SET UP LOG FILE */
			/* *************** */
			File logFile = new File(inputFolder + File.separator + "logs");
			logFile.mkdirs();
			logFile = new File(logFile + File.separator + name + ".log");
			MyPrintStream log = new MyPrintStream(logFile);
			
			/* ************** */
			/* SET PARAMETERS */
			/* ************** */
			DistortTetrahedra distort = new DistortTetrahedra();
			Thread_DistortTetrahedra distort_thread = new Thread_DistortTetrahedra(distort);
			distort.setInputFile(f);
			distort.setTimeStep(timeStep);
			distort.setLattice(lattice);
			distort.addIntraMolecPotential(new PotentialLookup(lookupPrecision, hookeCBrIntramolec));
			distort.addIntraMolecPotential(new PotentialLookup(lookupPrecision, hookeBrBrIntramolec));
			distort.addInterMolecPotential(new PotentialLookup(lookupPrecision, ljp_intermolec));
			distort.setIntermolecularMaxDist(intermolecularMaxDist);
			distort.setEndingBecauseOfEnergy(isEndingBecauseOfEnergy);
			distort.setEndDeltaE(endDeltaE);
			distort.setEndingBecauseOfLargeMovement(isEndingBecauseOfLargeMovement);
			distort.setLargeEndMovement(largeEndMovement);
			distort.setEndingBecauseOfSmallMovement(isEndingBecauseOfSmallMovement);
			distort.setSmallEndMovement(smallEndMovement);
			distort.setEndingBecauseOfTimeSteps(isEndingBecauseOfTimeSteps);
			distort.setEndTimeStep(endTimeStep);
			distort.setLog(log);
			distort_thread.setCheckEnergyEvery(stepsToCheckEnergy);
			
			/* ********************** */
			/* SET XYZ & MOVIE OUTPUT */
			/* ********************** */
			if(outputXYZs || outputMovie) {
				distort.setOutputtingXYZs(outputXYZs);
				distort.setOutputtingMovie(outputMovie);
				distort.setOutputXYZsEveryThisManyTimeSteps(outputXYZsEveryThisManyTimeSteps);
				distort.setOutputMovieXYZsEveryThisManyTimeSteps(outputMovieXYZsEveryThisManyTimeSteps);
				outputFolder = new File(f.getParentFile() + File.separator + "xyzOut");
				outputFolder.mkdirs();
				distort.setOutputFolder(outputFolder);
				distort.setOutputFilePrefix(name);
			}
			/* ******************** */
			/* SUBMIT THREAD TO RUN */
			/* ******************** */
			futures.add(es.submit(distort_thread));
		}
		boolean keepAlive = false;
		do {
			keepAlive = false;
			for(Future<?> future : futures)
				if(!future.isDone())
					keepAlive = true;
			try {
				Thread.sleep(10000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		} while(keepAlive);
		es.shutdown();
	}
}
