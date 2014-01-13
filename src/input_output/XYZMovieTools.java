package input_output;

import io.MyFileInputStream;
import io.MyPrintStream;

import java.io.File;
import java.util.Scanner;

public class XYZMovieTools {

	public static void splitMovieFile(File movieFile, int intoHowMany) {
		MyFileInputStream mfis = new MyFileInputStream(movieFile);
		Scanner s = mfis.getScanner();
		
		int numTimeSteps = 0;
		String[] line = null;
		while(s.hasNextLine()) {
			line = s.nextLine().split("\t");
			if(line.length == 1 && line[0].length() > 1)
				System.out.println("Time step: " + numTimeSteps++ + " found.");
		}
		
		System.out.println("Number of time steps in file: " + numTimeSteps);
		int stepsPerFile = numTimeSteps / intoHowMany;
		int numInLastFile = numTimeSteps % intoHowMany;
		
		mfis.close();
		mfis = new MyFileInputStream(movieFile);
		s = mfis.getScanner();
		
		MyPrintStream mps;
		String root = movieFile.getName();
		root = root.substring(0, root.lastIndexOf("."));
		String fileSuffix = root.substring(root.lastIndexOf("."));
		
		System.out.println("Input file suffix: " + fileSuffix);
		File outFile;
		for(int i = 0; i < intoHowMany; i++) {
			outFile = new File(movieFile.getParentFile().getAbsolutePath() + File.separator +
					root + "." + (i*stepsPerFile) + "-" + ((i+1)*stepsPerFile) + fileSuffix);
			mps = new MyPrintStream(outFile);
			String linesToPrint = "";
			String curLine;
			int numTimeStepsFound = 0;
			int numLinesToPrint = 0;
			while(s.hasNextLine() && numTimeStepsFound < stepsPerFile) {
				curLine = s.nextLine();
				if(curLine.split("\t").length == 1)
					numTimeStepsFound++;
				linesToPrint += curLine + "\n";
				if(++numLinesToPrint > 1000) {
					mps.print(linesToPrint);
					linesToPrint = "";
					numLinesToPrint = 0;
				}
			}
			mps.print(linesToPrint);
			mps.flush();
			mps.close();
		}

		/* PRINT THE LAST BIT TO A NEW FILE */
		outFile = new File(movieFile.getParentFile().getAbsolutePath() + File.separator +
				root + "." + (intoHowMany*stepsPerFile) + "-" + (numTimeSteps) + fileSuffix);
		mps = new MyPrintStream(outFile);
		String linesToPrint = "";
		String curLine;
		int numTimeStepsFound = 0;
		int numLinesToPrint = 0;
		while(s.hasNextLine() && numTimeStepsFound < stepsPerFile) {
			curLine = s.nextLine();
			linesToPrint += curLine + "\n";
			if(++numLinesToPrint > 1000) {
				mps.print(linesToPrint);
				linesToPrint = "";
				numLinesToPrint = 0;
			}
		}
		mps.print(linesToPrint);
		mps.flush();
		mps.close();
	}
	
	public static void main(String[] args) {
		File movieFile = new File("Z:\\Simulation\\Eric\\CBr4\\Dill-4\\89\\obj\\xyzOut\\movies\\1.mod.objmovie.xyz");
		int splitInto = 10;
		XYZMovieTools.splitMovieFile(movieFile, splitInto);
	}
}
