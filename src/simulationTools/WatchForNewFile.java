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
package simulationTools;

import input_output.LammpsTools;

import java.io.File;
import java.util.Observable;
import java.util.Observer;
import java.util.Vector;

import simulationTypes.RunSimulationThread;
import watchers.FolderWatcher;
import watchers.InvalidWatchSelectionException;

public class WatchForNewFile extends Observable implements Observer {
	private FolderWatcher fw;
	public final static String NEW_FILES_TO_COMPUTE = "New files for xray diffraction detected";
	public WatchForNewFile(File folderToWatch) {
		try {
			fw = new FolderWatcher(folderToWatch, FolderWatcher.WhatToWatchFor.NEW_FILE_IN_FOLDER);
		} catch (InvalidWatchSelectionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		fw.addObserver(this);
		fw.setMoveTo(new File(folderToWatch.getParentFile() + File.separator + "analyzed"));
		Thread t = new Thread(fw);
		t.start();
	}
	
	private File[] parse(File f) {
		LammpsTools.setNumUnitCells(15);
		return LammpsTools.LaampsDumpToXYZs(f);
	}
	private void calcDiffraction(File[] f) {
		setChanged();
		notifyObservers(new Object[] {NEW_FILES_TO_COMPUTE, f});
	}
	@Override
	public void update(Observable arg0, Object arg1) {
		if(arg0 instanceof FolderWatcher) {
			FolderWatcher fw = (FolderWatcher) arg0;
			if(arg1 instanceof Object[]) {
				Object[] obj = (Object[]) arg1;
				if(obj[0] instanceof FolderWatcher.WhatToWatchFor) {
					File moveTo, folder;
					File[] newFiles;
					Vector<File> actualNewFiles;
					switch((FolderWatcher.WhatToWatchFor) obj[0]) {
					case FOLDER_CREATION:
						System.out.println("The folder was created!");
						break;
					case NEW_FILE_IN_FOLDER:
						moveTo = (File) obj[1];
						folder = fw.getFolder();
						newFiles = folder.listFiles();
						actualNewFiles = new Vector<File>();
						for(File f : folder.listFiles()) {
							if(!f.isDirectory())
								actualNewFiles.add(f);
						}
						newFiles = new File[actualNewFiles.size()];
						newFiles = actualNewFiles.toArray(newFiles);
						String type = newFiles.length == 1 ? " file has" : " files have";
						System.out.println(newFiles.length + type + " appeared. ");
						if(moveTo != null)
							for(File file: newFiles) {
								File newFileName = new File(moveTo + File.separator + file.getParentFile().getName() + File.separator);
								newFileName.mkdirs();
								newFileName = new File(newFileName + File.separator + file.getName());
								if(file.renameTo(newFileName))
									System.out.println("Moving: " + file.getAbsolutePath() + "\n\tto:" + newFileName.getAbsolutePath());
								File[] parsed = null;
								if(newFileName.getName().contains(".xyz"))
									parsed = new File[] {newFileName};
								else
									parsed = parse(newFileName);
								calcDiffraction(parsed);
							}
						break;
					case NEW_FOLDER_IN_FOLDER:
						moveTo = (File) obj[1];
						folder = fw.getFolder();
						newFiles = folder.listFiles();
						actualNewFiles = new Vector<File>();
						for(File f : folder.listFiles()) {
							if(f.isDirectory())
								actualNewFiles.add(f);
						}
						newFiles = new File[actualNewFiles.size()];
						newFiles = actualNewFiles.toArray(newFiles);
						System.out.println(newFiles.length + " file(s) have appeared. ");
					}
				}
			}
		}
	}

	public static void main(String[] args) {
		File folderToWatch = new File("Z:\\Simulation\\Eric\\CBr4\\Dill-4\\88\\EDD_4-88n\\LAMMPS\\outputAtomicPositions");
		WatchForNewFile fwt = new WatchForNewFile(folderToWatch);
		
		
//		int step = 5;
//		RunSimulationThread simul = new RunSimulationThread();
//		
//		simul.setInputFolder(new File("Z:\\Simulation\\Eric\\CBr4\\Dill-4\\88\\EDD_4-88n\\LAMMPS\\dump\\parsedDumpToTimeSteps"));
//		simul.setObjFiles(simul.getInputFolder().listFiles());
//		simul.setOutputFolder(new File(simul.getInputFolder() + File.separator + "xrayCalc"));
//		simul.setqStep(1. / 15.);
//		int startModifyIdx = 0;
//		int finishModifyIdx = simul.getObjFiles().length;
//		int numSteps = (int) (finishModifyIdx / step);
//		int j = 0;
//		for(j = 0; j < numSteps; j++) {
//			simul.setStartModifyIdx(j*step);
//			simul.setFinishModifyIdx((j+1)*+step);
//			runDiffractionCalc(simul, projections, SimulationToRun.CALC_DIFFRACTION_XYZ);
//		}
//		simul.setStartModifyIdx(j);
//		simul.setFinishModifyIdx(finishModifyIdx);
//		RunSimulationThread.runDiffractionCalc(simul, projections, SimulationToRun.CALC_DIFFRACTION_XYZ);
//
//		simul.getEs_simul().shutdown();
//		simul.getEs_diffraction().shutdown();
	}
	
}
