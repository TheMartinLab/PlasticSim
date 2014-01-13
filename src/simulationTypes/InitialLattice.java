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
import java.util.Vector;

import defaultPackage.IdealTetrahedron;

public class InitialLattice implements Serializable, Lattices {

	private static final long serialVersionUID = 3343070542985579470L;
	private IdealTetrahedron[][][] lattice;
	private double a;
	private int numUnitCellsPerAxis;
	private long randomNumberSeed;
	private String initLatticeObjFileExtension;
	private String objFileName;
	
	public IdealTetrahedron[][][] getLattice() {
		return lattice;
	}
	public double getA() {
		return a;
	}
	public int getNumUnitCellsPerAxis() {
		return numUnitCellsPerAxis;
	}
	public long getRandomNumberSeed() {
		return randomNumberSeed;
	}
	public String getInitLatticeObjFileExtension() {
		return initLatticeObjFileExtension;
	}
	public void setLattice(IdealTetrahedron[][][] lattice) {
		this.lattice = lattice;
	}
	public void setA(double a) {
		this.a = a;
	}
	public void setNumUnitCellsPerAxis(int numUnitCellsPerAxis) {
		this.numUnitCellsPerAxis = numUnitCellsPerAxis;
	}
	public void setRandomNumberSeed(long randomNumberSeed) {
		this.randomNumberSeed = randomNumberSeed;
	}
	public void setInitLatticeObjFileExtension(String initLatticeObjFileExtension) {
		this.initLatticeObjFileExtension = initLatticeObjFileExtension;
	}
	@Override
	public String[] getParamsStringArray() {
		Vector<String> params = new Vector<String>();
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
		return params.toArray(new String[params.size()]);
	}
	public String getObjFileName() {
		return objFileName;
	}
	public void setObjFileName(String objFileName) {
		this.objFileName = objFileName;
	}
	
	
}
