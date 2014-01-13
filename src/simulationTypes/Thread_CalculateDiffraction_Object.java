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

import java.util.Date;
import java.util.Observable;

import defaultPackage.JAtom;
import defaultPackage.JPixel;
import defaultPackage.JavaToC;

public class Thread_CalculateDiffraction_Object extends Observable implements Runnable {

	private CalculateDiffraction calc;
	public final static String COMPLETE = "DIFFRACTION_COMPLETE";
	private MyPrintStream mpsLog;
	
	public Thread_CalculateDiffraction_Object(CalculateDiffraction calc, MyPrintStream mps) {
		this.calc = calc;
		mpsLog = mps;
		writeToLog("\tInitialized Thread_CalculateDiffraction for input file: " + calc.getShortObjectFileName());
	}
	
	private void writeToLog(String s) {
		mpsLog.println("Input file: " + calc.getShortObjectFileName() + " " + s);
	}
	public void run() {
		writeToLog("\tBeginning to calculate diffraction for: " + calc.getShortObjectFileName());
		Date startCalc = new Date();
		writeToLog("\t\tAbout to read in modified lattice: " + calc.getShortObjectFileName());
		JAtom[] lattice = null;
		try {
			calc.readInLattice();
		} catch (Exception e) {
			lattice = calc.readInXYZLattice();
			writeToLog("\t\tConverted molecules to atoms.");
			e.printStackTrace();
		} 
		calc.setqStep(1. / (double) calc.getLattice().getNumUnitCellsPerAxis()); 
		writeToLog("\t\tRead in modified lattice: " + calc.getShortObjectFileName());
		Date endCalc = new Date();
		writeToLog("\t\tTotal time to read lattice: " + (endCalc.getTime() - startCalc.getTime()) / 1000 + " s.");
		if(lattice == null) {
			lattice = calc.getFractionalAtoms();
			writeToLog("\t\tConverted molecules to atoms.");
		}
		JPixel[] pixels = calc.getPixels();
		int[] elemTypes = calc.getElemTypes();
		
		startCalc = new Date();
		JavaToC.calcDiffraction(lattice, pixels, elemTypes);
		endCalc = new Date();
		writeToLog("\t\tCalculated diffraction of " + lattice.length + " atoms onto " + pixels.length + 
				" pixels. Total time: " + (endCalc.getTime() - startCalc.getTime()) / 1000 + " s.");
		
		setChanged();
		notifyObservers(COMPLETE);
	}

	public CalculateDiffraction getCalc() {
		return calc;
	}

	public void setCalc(CalculateDiffraction calc) {
		this.calc = calc;
	}

	public MyPrintStream getMpsLog() {
		return mpsLog;
	}

	public void setMpsLog(MyPrintStream mpsLog) {
		this.mpsLog = mpsLog;
	}

}
