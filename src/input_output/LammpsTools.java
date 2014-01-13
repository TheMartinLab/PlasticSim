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
package input_output;

import io.MyFileInputStream;
import io.MyPrintStream;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;
import java.util.Vector;

import simulationTools.LatticeTools;

public class LammpsTools {

	public static void objectToLAMMPS(File inputRoot, File objRoot, File fileLAMMPSRoot) {
		fileLAMMPSRoot.mkdirs();
		
		File[] inputFiles = objRoot.listFiles();
		for(File inputFile : inputFiles) {
			System.out.println("Converting Java object to LAMMPS object: " + inputFile.getAbsolutePath());
 			String inputName = inputFile.getName();
			inputName = inputName.substring(0, inputName.lastIndexOf("."));
			File outputFile = new File(fileLAMMPSRoot + File.separator + "data." + inputName);
			File f = new File(outputFile + ".atomic");
//			LatticeTools.LatticeToLAMMPS_Atomic(f, inputFile);
			LatticeTools.LatticeToLAMMPS_Molecular(new File(outputFile + ".molecular"), inputFile);
			System.out.println("File: " + f.getName() + " written to disk");
		}
	}

	public final static String LINE_SEPARATOR = "ITEM: TIMESTEP";
	public final static String LINE_NUM_ATOMS = "ITEM: NUMBER OF ATOMS";
	public final static String LINE_BOX_BOUNDS = "ITEM: BOX BOUNDS pp pp pp";
	public final static String LINE_ATOM_START = "ITEM: ATOMS id type xs ys zs";
	public final static int ATOM_START_LINE = 10;
	public final static String DUMP_ROOT = ".atom";
	public static int numUnitCells;
	/**
	 * 
	 * @param fileLAMMPSRoot if directory: parse directory for files containing the .atom extension and then
	 * <br> split all those files into individual .xyz time steps
	 * <br> if file: split into individual .xyz time steps
	 * @return
	 */
	public static File[] LaampsDumpToXYZs(File fileLAMMPSRoot) {
		File[] files;
		if(fileLAMMPSRoot.isDirectory()) {
			files = fileLAMMPSRoot.listFiles();
			Vector<File> dumpFiles = new Vector<File>();
			for(File f : files) {
				if(f.getName().contains(DUMP_ROOT))
					dumpFiles.add(f);
			}
			return LaampsDumpToXYZs(dumpFiles.toArray(new File[dumpFiles.size()]));
		} else {
			return LaampsDumpToXYZs(new File[] {fileLAMMPSRoot});
		}
	}
	public static File[] LaampsDumpToXYZs(File[] files) {
		Vector<File> xyzFiles = new Vector<File>();
		MyFileInputStream mfis = null;
		for(File f : files) {
			mfis = new MyFileInputStream(f);
			Scanner s = mfis.getScanner();
			while(s.hasNextLine()) {
				String line = s.nextLine();
				if(line.compareTo(LINE_SEPARATOR) == 0) {
					int timeStep = Integer.valueOf(s.nextLine());
					xyzFiles.add(writeXYZ(s, timeStep, f));
				}
			}
		}
		mfis.close();
		return xyzFiles.toArray(new File[xyzFiles.size()]);
	}
	public static File writeXYZ(Scanner s, int timeStep, File curFile) {
		File dir = new File(curFile.getParentFile() + File.separator + 
				"parsedDumpToTimeSteps");
		dir.mkdirs();
		File f = new File(dir + File.separator + curFile.getName() + "." + timeStep);
		try {
			f.createNewFile();
		} catch (IOException e) {
			e.printStackTrace();
		}
		MyPrintStream mps = new MyPrintStream(f);
		int lineIdx = 0;
		int totalLineIdx = 0;
		int printPerLine = 1000;
		int totalAtoms = 100000;
		String lines = "";
		while(s.hasNextLine() && totalLineIdx != totalAtoms) {
			String line = s.nextLine();
			if(line.compareTo(LINE_NUM_ATOMS) == 0) {
				totalAtoms = Integer.valueOf(s.nextLine());
				mps.println(totalAtoms);
				mps.println();
			} else if(line.compareTo(LINE_BOX_BOUNDS) == 0) {
				s.nextLine();
				s.nextLine();
				s.nextLine();
			} else if(line.compareTo(LINE_ATOM_START) == 0) {
				while(s.hasNextLine() && totalLineIdx != totalAtoms) {
					line = s.nextLine();
					String[] split = line.split(" ");
					int id = Integer.valueOf(split[0]);
					int type = Integer.valueOf(split[1]);
					double x = Double.valueOf(split[2]);
					double y = Double.valueOf(split[3]);
					double z = Double.valueOf(split[4]);
					switch(type) {
					case 1:
						type = 6;
						break;
					case 2:
						type = 35;
						break;
					}
					lines += type + "\t" + x*numUnitCells + "\t" + y*numUnitCells + "\t" + z*numUnitCells + "\n";
					if(++lineIdx % printPerLine == 0) {					
						mps.print(lines);
						lineIdx = 0;
						lines = "";
					}
					totalLineIdx++;
				}
			}
		}
		mps.print(lines);
		mps.flush();
		mps.close();
		return f;
	}
	public static void main(String[] args) {
		File inputRoot = new File("Z:\\Simulation\\Eric\\CBr4\\Dill-4\\89");
		File objRoot = new File(inputRoot + File.separator + "obj\\analyzed obj");
		File fileLAMMPSRoot = new File(inputRoot + File.separator + "LAMMPS\\inputXYZ");
		boolean mkdir = fileLAMMPSRoot.mkdirs();
		if(args.length > 0) {
			fileLAMMPSRoot = new File(args[0]);
		}
		numUnitCells = 15;
		objectToLAMMPS(inputRoot, objRoot, fileLAMMPSRoot);
		
//		LaampsDumpToXYZs(fileLAMMPSRoot);
		
	}
	public static int getNumUnitCells() {
		return numUnitCells;
	}
	public static void setNumUnitCells(int numUnitCells) {
		LammpsTools.numUnitCells = numUnitCells;
	}
	
}
