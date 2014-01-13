package input_output;

import io.MyFileInputStream;
import io.MyPrintStream;
import io.StringConverter;

import java.io.File;
import java.util.Scanner;
import java.util.Vector;

public class SumXrayFilesInFolder {

	private File inFolder;
	private File outFolder;
	public final static String XRAY_DATA_START = "xyI_column";
	private boolean oldFormat = false;
	public SumXrayFilesInFolder() {
		
	}
	
	private double[][] getData(Scanner s) {
		Vector<double[]> data = new Vector<double[]>();
		
		String[] splitLine;
		double x, y, I;
		while(s.hasNextLine()) {
			splitLine = s.nextLine().split("\t");
			x = Double.valueOf(splitLine[0]);
			y = Double.valueOf(splitLine[1]);
			I = Double.valueOf(splitLine[2]);
			data.add(new double[] {x, y, I});
		}
		
		return data.toArray(new double[data.size()][]);
	}
	public void run() {
		/* INPUT & OUTPUT */
		File outFile,
			folders[] = inFolder.listFiles(),
			files[];
		MyPrintStream mps;
		MyFileInputStream mfis;
		Scanner s;
		String line;
		String[] splitLine;
		
		/* DATA VARIABLES */
		double x, y, I;	
		int numLines;
		
		boolean foldersIsFile = false;
		for(File folder : folders) {
			System.out.println("Reading from folder: " + folder.getAbsolutePath());
			outFile = new File(outFolder + File.separator + folder.getParentFile().getName() + "--" + folder.getName() + ".xray");
			mps = new MyPrintStream(outFile);
			files = folder.listFiles();
			if(folder.isFile()) {
				files = folders;
				foldersIsFile = true;
			}

			double[][] summedXrayData = null;
			Vector<String> header = new Vector<String>();
			Vector<String> fileNames = new Vector<String>(); 
			for(File file : files) {
				System.out.println("\tParsing file: " + file.getAbsolutePath());
				fileNames.add(file.getAbsolutePath());
				header.add("File: " + file.getAbsolutePath());
				mfis = new MyFileInputStream(file);
				s = mfis.getScanner();
				if(!oldFormat) {
					while(s.hasNextLine()) {
						line = s.nextLine();
						if(line.compareTo(XRAY_DATA_START) == 0)
							break;
						else
							header.add(line);						
					}
				} else 
					s.nextLine(); // skip the line that says how many were summed
				
				
				s.nextLine(); // skip the blank line
				double[][] fileData = getData(s);
				if(summedXrayData == null) { 
					summedXrayData = new double[fileData.length][3];
					for(int i = 0; i < fileData.length; i++) {
						summedXrayData[i][0] = fileData[i][0];
						summedXrayData[i][1] = fileData[i][1];
						summedXrayData[i][2] = fileData[i][2];
					}
				} else {
					for(int i = 0; i < fileData.length; i++) {
						summedXrayData[i][2] += fileData[i][2];
					}
				}
				
				s.close();
				mfis.close();
			}
			
			mps.printCommentLine(" Summed file names: ");
			for(String fileName : fileNames)
				mps.printCommentLine("\t" + fileName);
			
			line = "";
			numLines = 0;
			// TODO output summed xray data for each file
			mps.println("Summed Xray Data:\n");
			for(double[] d : summedXrayData) {
				line += StringConverter.arrayToTabString(d) + "\n";
				if(++numLines % 1000 == 0) {
					mps.print(line);
					numLines = 0;
					line = "";
				}
			}
			mps.print(line);
			mps.printCommentLine();
			mps.printCommentLine(" File headers: ");
			mps.printCommentLine();
			for(String headerLine : header) {
				mps.printCommentLine(" " + headerLine);
			}
			mps.flush();
			mps.close();
			if(foldersIsFile)
				break;
		}
			
	}
	public static void main(String[] args) {
		String root = "D:\\Documents referenced in lab notebooks\\Dill-4\\135\\diffraction";
		root = "Z:\\Simulation\\Eric\\CBr4\\Dill-4\\145\\analyzed\\diffraction";
		root = "D:\\Documents referenced in lab notebooks\\Dill-4\\146\\analyzed\\diffraction";
		String[] projections = new String[] {"001", "011", "111" };
		for(String proj : projections) {
			File inFolder = new File(root + File.separator + proj);
	//		inFolder = new File("D:\\$research\\current\\eclipse projects\\old_JNI_files\\diffraction\\old version of calculation with new version of control");
			File outFolder = new File(inFolder.getParent() + File.separator + "summed");
			outFolder.mkdirs();
			boolean oldFormat = false;
			
			SumXrayFilesInFolder sumThem = new SumXrayFilesInFolder();
			sumThem.oldFormat = oldFormat;
			
			sumThem.inFolder = inFolder;
			sumThem.outFolder = outFolder;
			
			sumThem.run();
		}
	}
}
