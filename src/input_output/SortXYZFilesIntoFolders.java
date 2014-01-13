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

import indexing.AlphabeticIndexingSystem;

import java.io.File;
import java.util.Vector;

import formatting.FormatDouble;

public class SortXYZFilesIntoFolders {

	private void byIndex(File root, String startIndex, String lastIndex) {
		AlphabeticIndexingSystem ais = new AlphabeticIndexingSystem(startIndex);
		File folder;
		while(ais.getName().compareTo(lastIndex) != 0) {
			folder = new File(root + File.separator + ais.getName() + File.separator + "xyz");
			if(!folder.exists())
				folder.mkdirs();
			File[] files = root.listFiles();
			for(File f : files)
				if(f.getName().startsWith(ais.getName()))
					f.renameTo(new File(folder + File.separator + f.getName()));
			
			ais.update();
		}
	}
	private void byProjection(File root, String[] projections) {
		for(String proj : projections) {
			Vector<File> matched = new Vector<File>();
			File[] files = root.listFiles();
			for(File f : files) {
				if(f.isDirectory())
					continue;
				if(f.getName().contains(proj)) {
					matched.add(f);
				}
			}
			File outFolder = new File(root + File.separator + proj);
			if(!outFolder.exists())
				outFolder.mkdirs();
			for(File f : matched) {
				f.renameTo(new File(outFolder + File.separator + f.getName()));
			}
		}
	}
	private void byTimeStep(File root) {
		double step = 1;
		double start = 1;
		double stop = 50;
		File[] files;
		String name;
		File folder;
		FormatDouble fd = new FormatDouble(FormatDouble.DecimalPlaces.TWO);
		String timeStep;	
		for(; start <= stop; start += step) {
			timeStep = "" + Math.rint(start * 100.)/100.;
			folder = new File(root + File.separator + "t=" + timeStep);
			if(!folder.exists())
				folder.mkdirs();
			files = root.listFiles();
			name = "." + timeStep + ".xyz";
			for(File f : files) {
				if(f.getName().contains(name) && !f.isDirectory())
					f.renameTo(new File(folder + File.separator + f.getName()));
					
			}
			
		}
	}
	
	public static void main(String[] args) {
		File root = new File("D:\\Documents referenced in lab notebooks\\Dill-4\\128\\analyzed\\EDD_146-a (5 shell) diffraction");
		root = new File("Z:\\Simulation\\Eric\\CBr4\\Dill-4\\146\\analyzed\\EDD_146-a (5 shell) diffraction");
		SortXYZFilesIntoFolders sort = new SortXYZFilesIntoFolders();
		String[] projections = new String[] {"001", "011", "111"};
		
//		sort.byTimeStep(root);
		String startIndex = "a";
		String endIndex = "d";
		sort.byIndex(root, startIndex, endIndex);
//		sort.byProjection(root, projections);
//		for(String proj : projections) {
//			File tmpRoot = new File(root + File.separator + proj);
//			sort.byIndex(tmpRoot, startIndex, endIndex);
//		}
		
		
		
	}
}
