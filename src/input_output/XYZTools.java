package input_output;

import defaultPackage.JVector;
import io.MyFileInputStream;

import java.io.File;
import java.util.Scanner;
import java.util.Vector;

import defaultPackage.JAtom;

public class XYZTools {
	
	public static JAtom[] readCrystalCoordinates(File inFile) {
		MyFileInputStream mfis;
		Scanner s;
		String splitLine[], line;
		int numAtoms = 0, Z;
		double x, y, z;
		Vector<JAtom> atoms;
		JVector pos;
		int numLines;
		
		mfis = new MyFileInputStream(inFile);
		s = mfis.getScanner();
		atoms = new Vector<JAtom>();
		while(s.hasNextLine()) {
			line = s.nextLine();
			splitLine = line.split("\t");
			
			if(splitLine.length == 1 && splitLine[0].length() > 0)
				numAtoms = Integer.valueOf(splitLine[0]);
			else if(splitLine.length == 4) {
				Z = Integer.valueOf(splitLine[0]);
				x = Double.valueOf(splitLine[1]);
				y = Double.valueOf(splitLine[2]);
				z = Double.valueOf(splitLine[3]);
				pos = new JVector(x, y, z);
				atoms.add(new JAtom(Z, pos));
			}
		}
		s.close();
		mfis.close();
		return atoms.toArray(new JAtom[atoms.size()]);
	}
	public static JAtom[] readXYZCoordinates(File inFile, double a) {
		MyFileInputStream mfis;
		Scanner s;
		String splitLine[], line;
		int numAtoms = 0, Z;
		double x, y, z;
		Vector<JAtom> atoms;
		JVector pos;
		int numLines;
		
		mfis = new MyFileInputStream(inFile);
		s = mfis.getScanner();
		atoms = new Vector<JAtom>();
		while(s.hasNextLine()) {
			line = s.nextLine();
			splitLine = line.split("\t");
			
			if(splitLine.length == 1 && splitLine[0].length() > 0)
				numAtoms = Integer.valueOf(splitLine[0]);
			else if(splitLine.length == 4) {
				Z = Integer.valueOf(splitLine[0]);
				x = Double.valueOf(splitLine[1]);
				y = Double.valueOf(splitLine[2]);
				z = Double.valueOf(splitLine[3]);
				pos = new JVector(x, y, z);
				pos.multiply(1./a);
				atoms.add(new JAtom(Z, pos));
			}
		}
		s.close();
		mfis.close();
		return atoms.toArray(new JAtom[atoms.size()]);
	}
}
