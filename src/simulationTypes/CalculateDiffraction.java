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

import io.MyFileInputStream;
import io.StringConverter;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;
import java.util.Vector;

import chemistry.JAtomTools;
import simulationTools.LatticeTools;
import simulationTypes.RunSimulationThread.DiffractionPrintStyle;
import simulationTypes.RunSimulationThread.Projections;
import defaultPackage.JPixel;
import defaultPackage.JVector;
import defaultPackage.JAtom;

public class CalculateDiffraction {
	private Lattices lattice;
	private Projections diffCalcType;
	private DiffractionPrintStyle diffractionPrintStyle;
	private JVector[][] axes;
	private double qMaxX;
	private double qMaxY;
	private double qStep;
	private JPixel[] pixels;
	private int[] elemTypes;
	private double wavelength;
	private double[][] xyI;
	private String objectFileName, shortObjectFileName;
	private File diffractionOutputFolder;
	private String index;
	private boolean areAtomsInCartesianCoordinates;
	
	public CalculateDiffraction(String objectFileName) {
		this.objectFileName = objectFileName;
		shortObjectFileName = objectFileName.substring(objectFileName.lastIndexOf(File.separator));
	}
	
	public void readInLattice() throws ClassNotFoundException, IOException {
		lattice = LatticeTools.readInLattice(new File(objectFileName));
	}
	
	public JAtom[] readInXYZLattice() {
		MyFileInputStream mfis = new MyFileInputStream(new File(objectFileName));
		Scanner s = mfis.getScanner();
		double maxCoordVal = 1.2 * lattice.getNumUnitCellsPerAxis();
		s.nextLine();
		Vector<JAtom> atoms = new Vector<JAtom>();
		JVector pos;
		while(s.hasNextLine()) {
			String[] splitLine = s.nextLine().split("\t");
			if(splitLine.length == 4) {
				int Z = Integer.valueOf(splitLine[0]);
				double x = Double.valueOf(splitLine[1]);
				double y = Double.valueOf(splitLine[2]);
				double z = Double.valueOf(splitLine[3]);
				pos = new JVector(x, y, z);
				if(areAtomsInCartesianCoordinates)
					pos.multiply(1./lattice.getA());
				
				atoms.add(new JAtom(Z, pos));
				
				if(!areAtomsInCartesianCoordinates && (x > maxCoordVal || y > maxCoordVal || z > maxCoordVal)) {
					areAtomsInCartesianCoordinates = true;
					atoms = null;
					mfis.close();
					return readInXYZLattice();
				}
			}
		}
		
		mfis.close();
		
		return atoms.toArray(new JAtom[atoms.size()]);
	}

	public void pixelsToDoubleArray() {
		for(int i = 0; i < pixels.length/xyI.length; i++) {
			int offset = i*xyI.length;
			for(int j = 0; j < xyI.length; j++) {
				xyI[j][2] += pixels[j + offset].getI(); 
			}
		}
	}
	
	/* ************** */
	/* STRING METHODS */
	/* ************** */
	public String[] getDiffractionParamsStringArray() {
		Vector<String> params = new Vector<String>();
		params.add("\tDiffraction parameters: ");
		params.add("\t\tDiffraction calculation type: " + diffCalcType.name());
		params.add("\t\tDiffraction axes: ");
		int idx = 0;
		for(JVector[] vec : axes)
			params.add("\t\t\t" + ++idx + ": " + StringConverter.arrayToTabString(vec));
		
		params.add("\t\tMaximum q values (x, y): " + qMaxX + "  " + qMaxY);
		params.add("\t\tDelta Q: " + qStep);
		params.add("\t\tElement Types: ");
		for(int Z : elemTypes) 
			params.add("\t\t\t" + JAtomTools.getAbbreviation(Z) + " " + JAtomTools.getName(Z));
		params.add("\t\tIncident wavelength: " + wavelength);
		params.add("\t\tNumber of pixels calculated: " + pixels.length);
		
		return params.toArray(new String[params.size()]);
	}
	public String[] getModifiedLatticeParamsArray() {
		return lattice.getParamsStringArray();
	}
	public JAtom[] getFractionalAtoms() {
		return LatticeTools.molToFractionalAtoms(lattice.getLattice(), lattice.getA());
	}
	public void nullify() {
		lattice = null;
		diffCalcType = null;
		diffractionPrintStyle = null;
		axes = null;
		pixels = null;
		elemTypes = null;
		xyI = null;
		objectFileName = null;
		shortObjectFileName = null;
		diffractionOutputFolder = null;
		index = null;
	}
	/* ********************************************* */
	/* GETTERS AND SETTERS FOR ALL PRIVATE VARIABLES */
	/* ********************************************* */
	public Projections getDiffCalcType() {
		return diffCalcType;
	}
	public void setDiffCalcType(Projections diffCalcType) {
		this.diffCalcType = diffCalcType;
	}
	public DiffractionPrintStyle getDiffractionPrintStyle() {
		return diffractionPrintStyle;
	}
	public void setDiffractionPrintStyle(DiffractionPrintStyle diffractionPrintStyle) {
		this.diffractionPrintStyle = diffractionPrintStyle;
	}
	public JVector[][] getAxes() {
		return axes;
	}
	public void setAxes(JVector[][] axes) {
		this.axes = axes;
	}
	public double getqMaxX() {
		return qMaxX;
	}
	public void setqMaxX(double qMaxX) {
		this.qMaxX = qMaxX;
	}
	public double getqMaxY() {
		return qMaxY;
	}
	public void setqMaxY(double qMaxY) {
		this.qMaxY = qMaxY;
	}
	public double getqStep() {
		return qStep;
	}
	public void setqStep(double qStep) {
		this.qStep = qStep;
	}
	public JPixel[] getPixels() {
		return pixels;
	}
	public void setPixels(JPixel[] pixels) {
		this.pixels = pixels;
	}
	public int[] getElemTypes() {
		return elemTypes;
	}
	public void setElemTypes(int[] elemTypes) {
		this.elemTypes = elemTypes;
	}
	public double getWavelength() {
		return wavelength;
	}
	public void setWavelength(double wavelength) {
		this.wavelength = wavelength;
	}
	public double[][] getXyI() {
		return xyI;
	}
	public void setXyI(double[][] xyI) {
		this.xyI = xyI;
	}
	public String getObjectFileName() {
		return objectFileName;
	}
	public void setObjectFileName(String objectFileName) {
		this.objectFileName = objectFileName;
	}
	public Lattices getLattice() {
		return lattice;
	}
	public void setLattice(Lattices lattice) {
		this.lattice = lattice;
	}

	public String getShortObjectFileName() {
		return shortObjectFileName;
	}

	public void setShortObjectFileName(String shortObjectFileName) {
		this.shortObjectFileName = shortObjectFileName;
	}

	public File getDiffractionOutputFolder() {
		return diffractionOutputFolder;
	}

	public void setDiffractionOutputFolder(File diffractionOutputFolder) {
		this.diffractionOutputFolder = diffractionOutputFolder;
	}

	public String getIndex() {
		return index;
	}

	public void setIndex(String index) {
		this.index = index;
	}

	public boolean areAtomsInCartesianCoordinates() {
		return areAtomsInCartesianCoordinates;
	}

	public void setareAtomsInCartesianCoordinates(
			boolean areAtomsInCartesianCoordinates) {
		this.areAtomsInCartesianCoordinates = areAtomsInCartesianCoordinates;
	}
}
