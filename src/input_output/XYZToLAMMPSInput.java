package input_output;

import geometry.JVector;
import io.MyFileInputStream;
import io.MyPrintStream;

import java.io.File;
import java.util.Scanner;
import java.util.Vector;

import chemistry.JAtom;
import chemistry.JAtomTools;

public class XYZToLAMMPSInput {

	private File inFolder, outFolder;
	private int[] ZToIdx;
	
	public void run() {
	
		ZToIdx = new int[100];
		
		MyFileInputStream mfis;
		Scanner s;
		
		MyPrintStream mps;
		String fileRoot;
		
		String line, splitLine[];
		
		int numAtoms, Z;
		double x, y, z;
		double xMin, xMax, yMin, yMax, zMin, zMax;
		
		Vector<JAtom> atoms;
		
		/* FIGURE OUT NUMBER OF ATOM TYPES */
		mfis = new MyFileInputStream(inFolder.listFiles()[0]);
		s = mfis.getScanner();
		numAtoms = Integer.valueOf(s.nextLine());
		s.nextLine();
		Vector<Integer> numZ = new Vector<Integer>();
		JVector pos;
		
		while(s.hasNextLine()) {
			splitLine = s.nextLine().split("\t");
			if(splitLine.length == 4) {
				Z = Integer.valueOf(splitLine[0]);
				x = Double.valueOf(splitLine[1]);
				y = Double.valueOf(splitLine[2]);
				z = Double.valueOf(splitLine[3]);
				if(!numZ.contains(Z)) {
					numZ.add(Z);
					ZToIdx[Z] = numZ.size();
				}
			}
		}
		
		s.close();
		mfis.close();
		String header = "";
		header += "atom_style atomic\n\n";
		header += numAtoms + " atoms\n\n";
		header += numZ.size() + " atom types\n\n";
		String masses = "";
		masses += "Masses\n\n";
		for(int i = 0; i < numZ.size(); i++)
			masses += (i+1) + " " + JAtomTools.getMass(numZ.get(i)) + "\n";

		
		
		
		
		for(File f : inFolder.listFiles()) {
			System.out.println("Reading: " + f.getAbsolutePath());
			mfis = new MyFileInputStream(f);
			s = mfis.getScanner();
			mps = new MyPrintStream(new File(outFolder + File.separator + f.getName() + ".data"));
			mps.print(header);
			atoms = new Vector<JAtom>();
			xMin = yMin = zMin = Double.MAX_VALUE;
			xMax = yMax = zMax = Double.MIN_VALUE;
			
			/* SKIP HEADER OF XYZ FILE */

			while(s.hasNextLine()) {
				splitLine = s.nextLine().split("\t");
				if(splitLine.length == 4) {
					Z = Integer.valueOf(splitLine[0]);
					Z = ZToIdx[Z];
					x = Double.valueOf(splitLine[1]);
					y = Double.valueOf(splitLine[2]);
					z = Double.valueOf(splitLine[3]);
					atoms.add(new JAtom(Z, new JVector(x, y, z)));
					
					if(xMin > x)
						xMin = x;
					if(yMin > y)
						yMin = y;
					if(zMin > z)
						zMin = z;
					
					if(xMax < x)
						xMax = x;
					if(yMax < y)
						yMax = y;
					if(zMax < z)
						zMax = z;
				}
			}
			s.close();
			mfis.close();
			
			mps.println(xMin + " " + xMax + " xlo xhi");
			mps.println(yMin + " " + yMax + " ylo yhi");
			mps.println(zMin + " " + zMax + " zlo zhi");
			mps.println();
			mps.println(masses);
			mps.println("Atoms");
			mps.println();
			int atomIdx = 0;
			line = "";
			for(JAtom atom : atoms) {
				pos = atom.getPosition();
				line += ++atomIdx + " " + atom.getZ() + " " + pos.i + " " + pos.j + " " + pos.k + "\n";
				if(atomIdx % 1000 == 0) {
					mps.print(line);
					line = "";
				}
			}
			mps.print(line);
			mps.flush();
			mps.close();
		}
		
		
	}
	
	public static void main(String[] args) {
		File inFolder = new File("Z:\\Simulation\\Eric\\CBr4\\Dill-4\\89\\obj\\analyzed\\xyzOut\\t=6.2\\in");
		File outFolder = new File(inFolder.getParent() + File.separator + "data.lammps");
		outFolder.mkdirs();
		
		XYZToLAMMPSInput xyz = new XYZToLAMMPSInput();
		xyz.inFolder = inFolder;
		xyz.outFolder = outFolder;

		xyz.run();
	}
}
