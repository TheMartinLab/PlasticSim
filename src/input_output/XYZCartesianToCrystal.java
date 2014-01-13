package input_output;

import geometry.JVector;
import io.MyFileInputStream;
import io.MyPrintStream;

import java.io.File;
import java.util.Scanner;
import java.util.Vector;

import chemistry.JAtom;

public class XYZCartesianToCrystal {

	private File inFolder, outFolder;
	private double a;
	public XYZCartesianToCrystal() {
		
	}
	public static JAtom[] XYZToCrystal(File inFile, double a) {
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
	public void run() {
		MyFileInputStream mfis;
		Scanner s;
		String splitLine[], line;
		MyPrintStream mps;
		
		int numAtoms = 0, Z;
		double x, y, z;
		Vector<JAtom> atoms;
		JVector pos;
		int numLines;
		
		for(File in : inFolder.listFiles()) {
			if(in.isDirectory())
				continue;
			System.out.println("Parsing file: " + in.getAbsolutePath());
			mfis = new MyFileInputStream(in);
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
			File out = new File(outFolder + File.separator + in.getName());
			mps = new MyPrintStream(out);

			System.out.println("Printing file: " + out.getAbsolutePath());
			mps.println(numAtoms + "\n");
			numLines = 0;
			line = "";
			JAtom atom;
			for(int i = 0; i < atoms.size(); i++) {
				atom = atoms.get(i);
				line += atom.getZ() + "\t" + atom.getPosition().toTabString() + "\n";
				if(++numLines % 1000 == 0) {
					mps.print(line);
					line = "";
					numLines = 0;
				}
			}
			mps.print(line);
			mps.flush();
			mps.close();
			
		}
	}
	public static void main(String[] args) {
		File inFolder = new File("Z:\\Simulation\\Eric\\CBr4\\Dill-4\\129\\simulOut\\cartesian xyz");
		File outFolder = new File(inFolder + File.separator + "xyzParsed");
		outFolder.mkdirs();
		double a = 1./8.82;
		
		XYZCartesianToCrystal xyz = new XYZCartesianToCrystal();
		
		xyz.inFolder = inFolder;
		xyz.outFolder = outFolder;
		xyz.a = a;
		
		xyz.run();
		
	}
}
