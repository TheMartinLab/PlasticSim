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

import java.io.Serializable;
import java.util.Iterator;
import java.util.Stack;
import java.util.Vector;

import chemistry.JAtomTools;
import defaultPackage.IdealTetrahedron;
import defaultPackage.LennardJonesPotential;

public class ModifiedLattice implements Serializable, Lattices {

	private static final long serialVersionUID = 3400674578948368035L;
	
	private IdealTetrahedron[][][] lattice;
	private double a;
	private int numUnitCellsPerAxis;
	private long randomNumberSeed;

	/* energy minimization strategy */
	private ModifyLatticeTypes modLattice;
	private String modifiedLatticeObjFileName;
	private String inputLatticeObjFileName;
	private Vector<String> previousIterationFileNames;

	private String outputLatticeObjFileName;
	private String outputLatticeXYZFileName;
	private double minimumIntermolecularDistance;
	private int threadIdx;
	/* interaction parameters */
	private double lookupPrecision;
	private LennardJonesPotential ljp;
	private int Z1;
	private int Z2;
	private double Z1Z2_distance;
	private double zeroPotentialDistance;
	private double wellDepth;
	private String energyUnits;
	
	/* random walk parameters */
	private int numRandomWalksPerMolecule;
	private int numStepsPerWalk;
	private double maximumTranslationalMotionPerMove;
	private double maximumRotationalAnglePerRotation;
	private double energyCutoff;
	private boolean useEnergyCutoff;
	private boolean disorderBetweenCycles;
	
	/* monte carlo orientational changes parameters */
	private int numTestsPerMolecule;
	private boolean maxNumChangesEnabled;
	private int maxNumOrientationalChanges;
	
	public ModifiedLattice() {
		previousIterationFileNames = new Stack<String>();
	}
	
	public void setCommonLatticeParams(Lattices lattice) {
		a = lattice.getA();
		this.lattice = lattice.getLattice();
		numUnitCellsPerAxis = lattice.getNumUnitCellsPerAxis();
		if(lattice.getObjFileName() != null)
			previousIterationFileNames.add(0, lattice.getObjFileName());
	}

	@Override
	public String[] getParamsStringArray() {
		Vector<String> params = new Vector<String>();
		
		/* GENERAL PARAMETERS */
		params.add("\tLattice parameter a: " + a);
		params.add("\tNumber of unit cells per axis: " + numUnitCellsPerAxis);
		int numMols = 0;
		for(int i = 0; i < lattice.length; i++)
			for(int j = 0; j < lattice[i].length; j++)
				for(int k = 0; k < lattice[i][j].length; k++)
					if(lattice[i][j][k] != null)
						numMols++;
		
		params.add("\tNumber of molecules in this lattice: " + numMols);
		params.add("\tRandom number seed: " + randomNumberSeed);
		params.add("\tModified lattice type: " + modLattice.name());
		params.add("\tPrevious Lattice input object file locations:");
		int idx = 0;
		for(String str : previousIterationFileNames) {
			params.add("\t\t" + (++idx) + ": " + str);
		}
		params.add("\tOutput lattice object file location: " + outputLatticeObjFileName);
		params.add("\tOutput lattice xyz file location: " + outputLatticeXYZFileName);
		params.add("\tCutoff intermolecular distance: " + minimumIntermolecularDistance + " (program will throw error if value is less than this)");
		params.add("\tThread index: " + threadIdx);
		
		/* LENNARD JONES INTERACTION PARAMETERS */
		params.add("\tLennard Jones Interaction Parameters: ");
		params.add("\t\tLookup precision: " + lookupPrecision + " (minimum distance precision");
		params.add("\t\tAtom type 1: " + JAtomTools.getName(Z1));
		params.add("\t\tAtom type 2: " + JAtomTools.getName(Z2));
		params.add("\t\tOptimum " + JAtomTools.getAbbreviation(Z1) + "-" + JAtomTools.getAbbreviation(Z2) +
				" distance: " + Z1Z2_distance);
		params.add("\t\tSigma parameter: " + zeroPotentialDistance + " (distance where LJ potential is zero");
		params.add("\t\tEpsilon parameter: " + wellDepth + " (Depth of LJ potential at optimum intermolecular distance");
		params.add("\t\tEnergy units: " + energyUnits);
		
		/* ENERGY MINIMIZATION PARAMETERS */
		params.add("\tGeneral energy minimization parameters: ");
		params.add("\t\tMaximum translational distance allowed per attempt: " + maximumTranslationalMotionPerMove);
		params.add("\t\tMaximum rotational angle allowed per attempt: " + maximumRotationalAnglePerRotation);
		params.add("\t\tEnergy cutoff: " + energyCutoff + " (Only move molecules that pairwise interactions greater than this value)");
		params.add("\t\tUse energy cutoff parameter? " + useEnergyCutoff);
		params.add("\t\tDisorder the lattice between cycles? " + disorderBetweenCycles);
		params.add("\tEnergy minimization strategy: " + modLattice.name());
		switch(modLattice) {
		case MONTE_CARLO_ORIENTATIONAL:
			params.add("\t\tNumber of monte carlo tests per molecule: " + numTestsPerMolecule);
			params.add("\t\tMaximum number of orientational changes per molecule: " + maxNumOrientationalChanges + 
					" (Energy minimization will exit after this number is exceeded.");
			params.add("\t\tMaximum number of orientational changes enabled? " + maxNumChangesEnabled);
			break;
		case RANDOM_WALK:
			params.add("\t\tNumber of random walks per molecule: " + numRandomWalksPerMolecule);
			params.add("\t\tNumber of steps per new random walk: " + numStepsPerWalk);
			break;
		}
		
		return params.toArray(new String[params.size()]);
	}

	public double getA() { return a; }
	public void setA(double a) { this.a = a; }
	public int getNumUnitCellsPerAxis() { return numUnitCellsPerAxis; }
	public void setNumUnitCellsPerAxis(int numUnitCellsPerAxis) { this.numUnitCellsPerAxis = numUnitCellsPerAxis; }
	public long getRandomNumberSeed() { return randomNumberSeed; }
	public void setRandomNumberSeed(long randomNumberSeed) { this.randomNumberSeed = randomNumberSeed; }
	public ModifyLatticeTypes getModLattice() { return modLattice; }
	public void setModLattice(ModifyLatticeTypes modLattice) { this.modLattice = modLattice; }
	public String getModifiedLatticeObjFileName() { return modifiedLatticeObjFileName; }
	public void setModifiedLatticeObjFileName(String modifiedLatticeObjFileName) { this.modifiedLatticeObjFileName = modifiedLatticeObjFileName; }
	public double getMinimumIntermolecularDistance() { return minimumIntermolecularDistance; }
	public void setMinimumIntermolecularDistance( double minimumIntermolecularDistance) { this.minimumIntermolecularDistance = minimumIntermolecularDistance; }
	public double getLookupPrecision() { return lookupPrecision; }
	public void setLookupPrecision(double lookupPrecision) { this.lookupPrecision = lookupPrecision; }
	public LennardJonesPotential getLjp() { return ljp; }
	public void setLjp(LennardJonesPotential ljp) { this.ljp = ljp; }
	public int getZ1() { return Z1; }
	public void setZ1(int z1) { Z1 = z1; }
	public int getZ2() { return Z2; }
	public void setZ2(int z2) { Z2 = z2; }
	public double getZ1Z2_distance() { return Z1Z2_distance; }
	public void setZ1Z2_distance(double z1z2_distance) { Z1Z2_distance = z1z2_distance; }
	public double getZeroPotentialDistance() { return zeroPotentialDistance; }
	public void setZeroPotentialDistance(double zeroPotentialDistance) { this.zeroPotentialDistance = zeroPotentialDistance; }
	public double getWellDepth() { return wellDepth; }
	public void setWellDepth(double wellDepth) { this.wellDepth = wellDepth; }
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
	public String getEnergyUnits() { return energyUnits; }
	public void setEnergyUnits(String energyUnits) { this.energyUnits = energyUnits; }
	@Override
	public IdealTetrahedron[][][] getLattice() { return lattice; }
	public void setLattice(IdealTetrahedron[][][] lattice) { this.lattice = lattice; }
	public String getInputLatticeObjFileName() { return inputLatticeObjFileName; }
	public void setInputLatticeObjFileName(String inputLatticeObjFileName) {
		previousIterationFileNames.add(0, inputLatticeObjFileName);
		this.inputLatticeObjFileName = inputLatticeObjFileName;
	}
	public String getOutputLatticeObjFileName() { return outputLatticeObjFileName; }
	public void setOutputLatticeObjFileName(String outputLatticeObjFileName) { this.outputLatticeObjFileName = outputLatticeObjFileName; }
	public String getOutputLatticeXYZFileName() { return outputLatticeXYZFileName; }
	public void setOutputLatticeXYZFileName(String outputLatticeXYZFileName) { this.outputLatticeXYZFileName = outputLatticeXYZFileName; }
	public int getThreadIdx() { return threadIdx; }
	public void setThreadIdx(int threadIdx) { this.threadIdx = threadIdx; }
	public Vector<String> getPreviousIterationFileNames() { return previousIterationFileNames; }
	public void setPreviousIterationFileNames(Vector<String> previousIterationFileNames) { this.previousIterationFileNames = previousIterationFileNames; }
	public void addPreviousIterationFileName(String fName) { previousIterationFileNames.add(0, fName); }
	@Override
	public String getObjFileName() { return getOutputLatticeObjFileName(); }
}
