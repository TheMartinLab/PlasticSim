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
package defaultPackage;

import java.io.File;

import javax.swing.JFileChooser;


public class JavaToC {

	public static native JPixel[] calcDiffraction(JAtom[] lattice, JPixel[] pixels, int[] elemTypes);
	
	public static native void cuInit();
	
	static {
		String prop = "java.library.path";
		String curPath = System.getProperty(prop);
		String newPath = "D:\\$research\\current\\eclipse projects\\JNI\\";
		if(!curPath.contains(newPath)) {
			System.out.println(prop + " before setting to: " + newPath + ": " + System.getProperty(prop));
			System.setProperty(prop, curPath + ";" + newPath);
			System.out.println(System.getProperty(prop));
			System.out.println(prop + " after setting to: " + newPath + ": " + System.getProperty(prop));
		}
		String fileToLoad = "Diffraction64";
		if(!new File(fileToLoad + ".dll").exists()) {
			File loc = new File(".");
			System.out.println("Current file location: " + loc.getAbsolutePath());
			JFileChooser chooser = new JFileChooser(loc);
			chooser.setMultiSelectionEnabled(false);
			chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			int returnVal = chooser.showOpenDialog(null);
			switch(returnVal) {
			case JFileChooser.APPROVE_OPTION:
				fileToLoad = chooser.getSelectedFile().getName();
				fileToLoad = fileToLoad.substring(0, fileToLoad.lastIndexOf("."));
				break;
			}
		}
		
		System.loadLibrary(fileToLoad);
		
		//cuInit();
	}
}
