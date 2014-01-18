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
import java.util.Vector;

import defaultPackage.Potential;
import defaultPackage.PotentialLookup;


public class DistortTetrahedra {

	/* ***************** */
	/* LATTICE VARIABLES */
	/* ***************** */
	private Lattices lattice;
	private double tolerance = 1e-4;
	private double intermolecularMaxDist;
	private double timeStep;
	
	/* ********** */
	/* POTENTIALS */
	/* ********** */
	private Vector<PotentialLookup> interMolecPotentials;
	private Vector<PotentialLookup> intraMolecPotentials;
	/* ************************ */
	/* INPUT & OUTPUT VARIABLES */
	/* ************************ */
	private File inputFile;
	private File outputFolder;
	private boolean isOutputtingXYZs;
	private boolean isOutputtingMovie;
	private String outputFilePrefix;
	private MyPrintStream log;
	private int outputXYZsEveryThisManyTimeSteps;
	private int outputMovieXYZsEveryThisManyTimeSteps;
	
	/* ************** */
	/* END CONDITIONS */
	/* ************** */
	private boolean isEndingBecauseOfEnergy = false;
	private boolean isEndingBecauseOfSmallMovement = false;
	private boolean isEndingBecauseOfLargeMovement = false;
	private boolean isEndingBecauseOfTimeSteps = false;
	private double endDeltaE;
	private double smallEndMovement;
	private double largeEndMovement;
	private int endTimeStep;
	
	/* ************** */
	/* MISC VARIABLES */
	/* ************** */
	
	/* ************** */
	/* CONSTRUCTOR(S) */
	/* ************** */
	public DistortTetrahedra() {
		
	}
	/* ***************** */
	/* FINALIZING OBJECT */
	/* ***************** */
	public void close() {
		log.flush();
		log.close();
		outputFolder = null;
		interMolecPotentials = null;
		intraMolecPotentials = null;
	}
	/* ************** */
	/* STRING METHODS */
	/* ************** */
	public String[] getParamsArray() {
		Vector<String> params = new Vector<String>();
		
		params.add("\tDistortTetrahedra parameters");
		params.add("\t\tTime Step: " + timeStep);
		for(PotentialLookup potLook : interMolecPotentials) {
			params.add("\t\tPotential: " + potLook.getType());
			params.add("\t\t\tZ1 Z1: " + potLook.getZ1() + " " + potLook.getZ2());
			params.add("\t\t\tdist potential" + potLook.getDistVal() + " " + potLook.getPotVal());
		}
		return params.toArray(new String[params.size()]);
	}
	/* ****** */
	/* OUTPUT */
	/* ****** */
	public void logParameters() {
		/* LATTICE VARIABLES */
		log.println("Lattice Variables");
		for(String line : lattice.getParamsStringArray())
			log.println(line);
		log.println();
		log.println("\tMax intermolecular distance before boundary conditions kick in: " + intermolecularMaxDist);
		log.println("\tSimulation time step: " + timeStep);

		/* POTENTIALS */
		log.println();
		log.println("Intermolecular Potential(s)");
		for(PotentialLookup potLook : interMolecPotentials)
			log.println("\t" + potLook.toString());
		log.println();
		log.println("Intramolecular Potential(s)");
		for(PotentialLookup potLook : intraMolecPotentials)
			log.println("\t" + potLook.toString());

		/* INPUT & OUTPUT VARIABLES */
		log.println();
		log.println("Input & Output Variables");
		log.println("\tInput folder: " + inputFile.getAbsolutePath());
		log.println("\tOutput folder: " + outputFolder.getAbsolutePath());
		log.println("\tLog file: " + log.getFile().getAbsolutePath());
		log.print("\tOutputting .xyz files: " + isOutputtingXYZs);
		if(isOutputtingXYZs)
			log.println("\tEvery " + outputXYZsEveryThisManyTimeSteps + " time steps.");
		else
			log.println();
		log.print("\tOutputting an .xyz movie file: " + isOutputtingMovie);
		if(isOutputtingXYZs)
			log.println("\tEvery " + outputMovieXYZsEveryThisManyTimeSteps + " time steps.");
		else log.println();
		
		/* SIMULATION END CONDITIONS */
		log.println();
		log.println("Simulation End Conditions");
		log.print("\tEnding because of deltaE: " + isEndingBecauseOfEnergy);
		if(isEndingBecauseOfEnergy)
			log.println(" Value: " + endDeltaE);
		else
			log.println();
		log.print("\tEnding because minimum movement target achieved: " + isEndingBecauseOfSmallMovement);
		if(isEndingBecauseOfSmallMovement)
			log.println(" Value: " + smallEndMovement);
		else
			log.println();
		log.print("\tEnding because maximum movement allowable exceeded: " + isEndingBecauseOfLargeMovement);
		if(isEndingBecauseOfLargeMovement)
			log.println(" Value: " + largeEndMovement);
		else
			log.println();
		log.print("\tEnding because of time step: " + isEndingBecauseOfTimeSteps);
		if(isEndingBecauseOfTimeSteps)
			log.println(" Value: " + endTimeStep);
		else
			log.println();
	}
	/* ********** */
	/* POTENTIALS */
	/* ********** */
	public void addInterMolecPotential(PotentialLookup potLook) { 
		if(interMolecPotentials == null)
			interMolecPotentials = new Vector<PotentialLookup>();
		
		interMolecPotentials.add(potLook);
	}
	public void addIntraMolecPotential(PotentialLookup potLook) { 
		if(intraMolecPotentials == null)
			intraMolecPotentials = new Vector<PotentialLookup>();
		
		intraMolecPotentials.add(potLook);
	}
	/* ***************** */
	/* GETTERS & SETTERS */
	/* ***************** */
	/* LATTICE */
	public Lattices getLattice() { return lattice; }
	public void setLattice(Lattices lattice) { this.lattice = lattice; }
	public double getTolerance() { return tolerance; }
	public void setTolerance(double tolerance) { this.tolerance = tolerance; }
	public double getIntermolecularMaxDist() { return intermolecularMaxDist; }
	public void setIntermolecularMaxDist(double intermolecularMaxDist) { this.intermolecularMaxDist = intermolecularMaxDist; }
	public double getTimeStep() { return timeStep; }
	public void setTimeStep(double timeStep) { this.timeStep = timeStep; }
	/* POTENTIALS */
	public Vector<PotentialLookup> getInterMolecPotentials() { return interMolecPotentials; }
	public void setInterMolecPotentials(Vector<PotentialLookup> interMolecPotentials) { this.interMolecPotentials = interMolecPotentials; }
	public Vector<PotentialLookup> getIntraMolecPotentials() { return intraMolecPotentials; }
	public void setIntraMolecPotentials(Vector<PotentialLookup> intraMolecPotentials) { this.intraMolecPotentials = intraMolecPotentials; }
	/* INPUT & OUTPUT */
	public File getInputFile() { return inputFile; }
	public void setInputFile(File inputFile) { this.inputFile = inputFile; }
	public File getOutputFolder() { return outputFolder; }
	public void setOutputFolder(File outputFolder) { this.outputFolder = outputFolder; }
	public boolean isOutputtingXYZs() { return isOutputtingXYZs; }
	public void setOutputtingXYZs(boolean isOutputtingXYZs) { this.isOutputtingXYZs = isOutputtingXYZs; }
	public boolean isOutputtingMovie() { return isOutputtingMovie; }
	public void setOutputtingMovie(boolean isOutputtingMovie) { this.isOutputtingMovie = isOutputtingMovie; }
	public String getOutputFilePrefix() { return outputFilePrefix; }
	public void setOutputFilePrefix(String outputFilePrefix) { this.outputFilePrefix = outputFilePrefix; }
	public MyPrintStream getLog() { return log; }
	public void setLog(MyPrintStream log) { this.log = log; }
	public int getOutputXYZsEveryThisManyTimeSteps() { return outputXYZsEveryThisManyTimeSteps; }
	public void setOutputXYZsEveryThisManyTimeSteps(int outputXYZsEveryThisManyTimeSteps) { this.outputXYZsEveryThisManyTimeSteps = outputXYZsEveryThisManyTimeSteps; }
	public int getOutputMovieXYZsEveryThisManyTimeSteps() { return outputMovieXYZsEveryThisManyTimeSteps; }
	public void setOutputMovieXYZsEveryThisManyTimeSteps(int outputMovieXYZsEveryThisManyTimeSteps) { this.outputMovieXYZsEveryThisManyTimeSteps = outputMovieXYZsEveryThisManyTimeSteps; }
	/* END CONDITIONS */
	public boolean isEndingBecauseOfEnergy() { return isEndingBecauseOfEnergy; }
	public void setEndingBecauseOfEnergy(boolean isEndingBecauseOfEnergy) { this.isEndingBecauseOfEnergy = isEndingBecauseOfEnergy; }
	public boolean isEndingBecauseOfSmallMovement() { return isEndingBecauseOfSmallMovement; }
	public void setEndingBecauseOfSmallMovement(boolean isEndingBecauseOfSmallMovement) { this.isEndingBecauseOfSmallMovement = isEndingBecauseOfSmallMovement; }
	public boolean isEndingBecauseOfLargeMovement() { return isEndingBecauseOfLargeMovement; }
	public void setEndingBecauseOfLargeMovement( boolean isEndingBecauseOfLargeMovement) { this.isEndingBecauseOfLargeMovement = isEndingBecauseOfLargeMovement; }
	public double getEndDeltaE() { return endDeltaE; }
	public void setEndDeltaE(double endDeltaE) { this.endDeltaE = endDeltaE; }
	public double getSmallEndMovement() { return smallEndMovement; }
	public void setSmallEndMovement(double smallEndMovement) { this.smallEndMovement = smallEndMovement; }
	public double getLargeEndMovement() { return largeEndMovement; }
	public void setLargeEndMovement(double largeEndMovement) { this.largeEndMovement = largeEndMovement; }
	public boolean isEndingBecauseOfTimeSteps() { return isEndingBecauseOfTimeSteps; }
	public void setEndingBecauseOfTimeSteps(boolean isEndingBecauseOfTimeSteps) { this.isEndingBecauseOfTimeSteps = isEndingBecauseOfTimeSteps; }
	public int getEndTimeStep() { return endTimeStep; }
	public void setEndTimeStep(int endTimeStep) { this.endTimeStep = endTimeStep; }
	/* MISC VARIABLES */
}
