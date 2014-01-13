package input_output;

import geometry.JVector;
import io.MyPrintStream;

import java.io.File;

import chemistry.JAtom;

public class CollapseSimulation {

	public static void main(String[] args) {
		String runIdx = "f";
		String runRoot = "mod";
		double a = 1;
		int startIdx = 0;
		int endIdx = 20;
		int numFiles = endIdx - startIdx + 1;
		JAtom[][] arr = new JAtom[numFiles][];
		File inFolder = new File("D:\\Documents referenced in lab notebooks\\Dill-4\\130\\simulOut\\xyz\\" +
				runIdx + "\\xyz\\");
		File outFile = new File(inFolder + File.separator + 
				+ startIdx + "-" + endIdx + ".collapsed.xyz");
		MyPrintStream mps = new MyPrintStream(outFile);
		
		int numAtoms = 0;
		for(int i = 0; i < numFiles; i++) {
			File inFile = new File(inFolder + File.separator + runIdx + "--" + i + "--" + runRoot + ".xyz");
			System.out.println("Reading file " + (i+1) + " of " + numFiles + ": " + inFile.getAbsolutePath());
			arr[i] = XYZCartesianToCrystal.XYZToCrystal(inFile, a);
			numAtoms += arr[i].length;
		}

		mps.println(numAtoms + "\n");
		
		for(int i = 0; i < arr.length; i++) {
			System.out.println("Printing atoms " + (i*arr[i].length) + " to " + ((i+1) * arr[i].length) + " of " + numAtoms);
			JVector pos;
			int idx = 0;
			String line = "";
			JVector trans = new JVector(2, 2, 2);
			for(JAtom atom : arr[i]) {
				pos = atom.getPosition();
				pos.add(trans);
				pos.i %= 1;
				pos.j %= 1;
				pos.k %= 1;
				atom.setPosition(pos);
				line += atom.toStringForXYZ() + "\n";
				if(idx++ > 1000) {
					idx = 0;
					mps.print(line);
					line = "";
				}
			}
			mps.print(line);
		}
		mps.close();
		
	}
}
