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
/**
 * @author Eric D. Dill eddill@ncsu.edu
 * @author James D. Martin jdmartin@ncsu.edu
 * Copyright © 2010-2013 North Carolina State University. All rights reserved
 */

import defaultPackage.JVector;

import input_output.XYZTools;
import io.MyPrintStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.JFileChooser;

class Diffraction {
	/* instance variables */
	private final static boolean inEclipse = false;
	private double qStep, wavelength;
	private JPixel[] p001, p110, p111;
	private JPixel[][] pAll;
	private double[][] xyI001, xyI110, xyI111;
	private double[][] qI;
	private int[] elemTypes;
	private JAtom[] lattice;
	private LennardJonesPotential ljp;
	private PotentialLookup potLook;
	HashMap<Double, JComplex> sfMap;
	private static FileOutputStream fos;
	private static PrintStream ps;
	/* native methods */
	private boolean b100 = false;
	private boolean b110 = false;
	private boolean b111 = true;
	private boolean bAll = false;
	private static boolean bit32 = false;
	private double qxMax = 10;
	//private BravaisLattice b;
	private OneDimensionalDiffraction oneD = new OneDimensionalDiffraction();
	int[] totalWalks;
	/* class methods */
	public static void printLatticeXYZ(JAtom[] lattice, String fName) {
		FileOutputStream fos = null;
		PrintStream ps = null;
		try {
			fos = new FileOutputStream(fName);
			ps = new PrintStream(fos);
		} catch(FileNotFoundException fnfe) {
			System.out.println("File not found: " + fName);
			return;
		}
		ps.println(lattice.length+"\n");
		String lines = "";
		for(int i = 0; i < lattice.length; i++) {
			if(i%1000 == 0) {
				ps.print(lines);
				lines = "";
			}
			lines += lattice[i].toStringForXYZ();
		}
		ps.println(lines);
		ps.close();
		try { fos.close(); }
		catch(IOException e) {}
	}
	private double[][] coalesce(double[][] xyI) {
		HashMap<Double, double[]> pix = new HashMap<Double, double[]>(1000);
		double len;
		double[] vals;
		DecimalFormat twoDForm = new DecimalFormat("#.####");
		for(int a = 0; a < pAll.length; a++) {
			for(int b = 0; b < pAll[a].length; b++) {
				len = Double.valueOf(twoDForm.format(pAll[a][b].getQ().length()));
				if(pix.containsKey(len)) {
					vals = pix.get(len);
					vals[0] += pAll[a][b].getI();
					vals[1]++;
					pix.put(len, vals);
				} else {
					pix.put(len, new double[] {pAll[a][b].getI(), 1});
				}
			}
		}
		int mapSize = pix.size();
		double[][] qINum = new double[mapSize][3]; 
		Iterator<Double> keySet = pix.keySet().iterator();
		Double key;
		int idx = 0;
		while(keySet.hasNext()) {
			key = keySet.next();
			vals = pix.get(key);
			qINum[idx][0] = key;
			qINum[idx][1] = vals[0];
			qINum[idx][2] = vals[1];
			idx++;
		}
		return qINum;
	}
	public static void printPixelsToFile(double[][] xyI, String fName, int numTimes) {
		File f = new File(fName);
		f.getParentFile().mkdirs();
		MyPrintStream mps = new MyPrintStream(f);
		mps.println("numTimes: " + numTimes);
		String lines = "";
		for(int i = 0; i < xyI.length; i++) {
			if(i%1000 == 0) {
				mps.print(lines);
				lines = "";
			}
			switch(xyI[0].length) {
			case 2: 
				lines += "\n" + xyI[i][0] + "\t" + xyI[i][1];
				break;
			case 3:
				lines += "\n" + xyI[i][0] + "\t" + xyI[i][1] + "\t" + xyI[i][2];
				break;
			}
		}
		mps.print(lines);
		mps.close();
	}
	/* instance methods */
	public void addToXYI() {
		int offset = 0;
		if(b100) {
			for(int i = 0; i < p001.length/xyI001.length; i++) {
				offset = i*xyI001.length;
				for(int j = 0; j < xyI001.length; j++) {
					xyI001[j][2] += p001[j + offset].getI(); 
				}
			}
		}
		if(b110) {
			for(int i = 0; i < p110.length/xyI110.length; i++) {
				offset = i*xyI110.length;
				for(int j = 0; j < xyI110.length; j++) {
					xyI110[j][2] += p110[j + offset].getI(); 
				}
			}
		}
		if(b111) {
			for(int i = 0; i < p111.length/xyI111.length; i++) {
				offset = i*xyI111.length;
				for(int j = 0; j < xyI111.length; j++) {
					xyI111[j][2] += p111[j + offset].getI(); 
				}
			}
		}
	}
	public void addToXYI(JPixel[] pixels, double[][] xyI) {
		int offset = 0;
		for(int i = 0; i < pixels.length/xyI.length; i++) {
			offset = i*xyI.length;
			for(int j = 0; j < xyI.length; j++) {
				xyI[j][2] += pixels[j + offset].getI(); 
			}
		}
	}
	public JAtom[] molToAtoms(Molecule[][][] lattice, JVector axes) {
		int numAtoms = 0;
		for(int i = 0; i < lattice.length; i++)
			for(int j = 0; j < lattice[i].length; j++)
				for(int k = 0; k < lattice[i][j].length; k++)
					if(lattice[i][j][k] != null)
						numAtoms += lattice[i][j][k].getAtoms().length;
		JAtom[] temp = new JAtom[numAtoms];
		JAtom[] mol;
		int atomIdx = 0;
		for(int i = 0; i < lattice.length; i++){
			for(int j = 0; j < lattice[i].length; j++){
				for(int k = 0; k < lattice[i][j].length; k++){
					if(lattice[i][j][k] != null) {
						mol = lattice[i][j][k].getAtoms();
						for(int a = 0; a < mol.length; a++) {
							temp[atomIdx] = (JAtom) mol[a].clone();
							//if(temp[atomIdx].Z == 6) { continue; }
							temp[atomIdx].getPosition().i = mol[a].getPosition().i / axes.i;
							temp[atomIdx].getPosition().j = mol[a].getPosition().j / axes.j;
							temp[atomIdx].getPosition().k = mol[a].getPosition().k / axes.k;
							atomIdx++;
						}
					}
				}
			}
		}
		return temp;
	}
	public void cartesianToFractional(JAtom[] lattice, double aConstant) {
		for(int i = 0; i < lattice.length; i++) {
			lattice[i].setPosition(JVector.multiply(lattice[i].getPosition(), 1./aConstant));
		}
	}
	enum DiffractionType {
		_2d_001,
		_2d_010,
		_2d_100,
		_2d_110,
		_2d_111,
		;
	}
	DiffractionType type = DiffractionType._2d_001;
	public void calcDiffraction(JAtom[] lattice, JPixel[] pixels, double[][] xyI, 
			boolean isCartesianCoordinates, String outFileName) {
		
		// set up parameters
		System.out.println("Number of atoms: " + lattice.length);
		
		// zero pixel intensity
		for(JPixel p : pixels)
			p.setI(0);
		for(int i = 0; i < xyI.length; i++)
			xyI[i][2] = 0;
		
		// compute diffraction
		JavaToC.calcDiffraction(lattice, pixels, elemTypes);
		addToXYI(pixels, xyI);
		printPixelsToFile(xyI, outFileName, 0);
	}
	public void calcDiffractionWSimul(String root) {
		String fName001 = root + "001.xray";
		String fName110 = root + "110.xray";
		String fName111 = root + "111.xray";
		
		if(b100) {
			System.out.println("\n\nDiffraction: 001\n\n");
			JavaToC.calcDiffraction(lattice, p001, elemTypes);
			//ewaldProjection(p001);
		}
		if(b110) {
			System.out.println("\n\nDiffraction: 110\n\n");
			JavaToC.calcDiffraction(lattice, p110, elemTypes);
		}
		if(b111) {
			System.out.println("\n\nDiffraction: 111\n\n");
			JavaToC.calcDiffraction(lattice, p111, elemTypes);
		}
		
		addToXYI();
		if(b100) { 
			printPixelsToFile(xyI001, fName001, 0); 
		}
		if(b110) { 
			printPixelsToFile(xyI110, fName110, 0); 
		}
		if(b111) { 
			printPixelsToFile(xyI111, fName111, 0); 
		}
		
		zeroXYI();
		zeroPixels();
	}
	private void zeroPixels() {
		if(p001 != null)
			for(int i = 0; i < p001.length; i++)
				p001[i].setI(0);
		
		if(p110 != null)
			for(int i = 0; i < p110.length; i++)
				p110[i].setI(0);
		
		if(p111 != null)
			for(int i = 0; i < p111.length; i++)
				p111[i].setI(0);
	}
	private void zeroXYI() {
		if(xyI001 != null)
			for(int i = 0; i < xyI001.length; i++)
				xyI001[i][2] = 0;
		
		if(xyI110 != null)
			for(int i = 0; i < xyI110.length; i++)
				xyI110[i][2] = 0;
		
		if(xyI111 != null)
			for(int i = 0; i < xyI111.length; i++)
				xyI111[i][2] = 0;
	}
	public int calcDiffractionWSimul(int numTimes, int numUnits, int walkLength, int numWalks, int monteCarloOrientations, int monteCarloShifts,
			int whichBox, double percentMono, double percentSecond, double percentFirst, double enerCutoff, double ljDepth, int mdCycles, 
			double distCutoff, int box, double maxTrans, double maxTorque) {
		
		boolean disorderBetweenCycles = false;
		double a = 8.8;
		System.out.println("numWalks: " + numWalks + "\tWalkLength: " + walkLength);
		Lattice l = new Lattice(numUnits, a);
		Simulate s = new Simulate(potLook, 0, a);
		//JVector hydrateAxes = new JVector(6.55, 6.6, 15.12);
		JVector hydrateAxes = new JVector(a, a, a);
		s.setTorque(.25);
		s.setTrans(.25);
		l.setS(s);
		l.distanceCutoff = 3.0;
		int walks = 0;
		try {
			l.readBoxes(new File("sorted_eight.lattice"), 8);
			l.readBoxes(new File("sorted_rotated_second_shells.lattice"), 9);
			l.readBoxes(new File("1.5_truncatedRotatedBoxes.xyz"), whichBox);
			for(int i = 1; i <= 4; i++) {
				//l.readBoxes(new File(i + "1.5_truncatedBoxes.xyz"), i);/*
				/*switch(i) {
				case 1: if(l.boxes1[0][0] == null) { l.readBoxes(new File(i + ".0_truncatedBoxes.xyz"), i); break;}
				case 2: if(l.boxes2[0][0] == null) { l.readBoxes(new File(i + ".0_truncatedBoxes.xyz"), i); break;}
				case 3: if(l.boxes3[0][0] == null) { l.readBoxes(new File(i + ".0_truncatedBoxes.xyz"), i); break;}
				case 4: if(l.boxes4[0][0] == null) { l.readBoxes(new File(i + ".0_truncatedBoxes.xyz"), i); break;}
				case 5: if(l.boxes5[0][0] == null) { l.readBoxes(new File(i + ".0_truncatedBoxes.xyz"), i); break;}
				case 6: if(l.boxes6[0][0] == null) { l.readBoxes(new File(i + ".0_truncatedBoxes.xyz"), i); break;}
				case 7: if(l.boxes7[0][0] == null) { l.readBoxes(new File(i + ".0_truncatedBoxes.xyz"), i); break;}
				default: if(l.boxes2[0][0] == null) { l.readBoxes(new File(2 + ".0_truncatedBoxes.xyz"), i); break;}
				}*/
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.out.println("file not found");
		}
		Clock c, c1;
		String fName001, fName110, fName111, fName1D;
		String[] fixed = {"xy", "xz", "yz"};
		double[] amt = {.001};
		int scatterer = 35;
		JPixel[] p;
		int[][] xyI;
		String proj;
		Molecule[][][] mols = null;
		String prefix = "torque_trans-" + maxTorque + "_" + maxTrans;
		if(disorderBetweenCycles) { prefix += "_disorderBetweenCycles"; }
		s.setMaxTorque(maxTorque);
		s.setMaxTrans(maxTrans);
		prefix += "refreshingCBr4Stuff";
		//prefix += "_Tetragonal_DispersivePotential_";
		//prefix += "_Tetragonal_LJPotential_";
		//prefix += "_MonoBox_RandomLocations_" + box + "_of_32_";
		//prefix += "_hyperbolic_";
		//prefix = "tryFirstCoordinate_PairWalk+monteCarloFill_" + percentMono + "_" + monteCarloOrientations + "_" + monteCarloShifts + "_" ;
		//prefix += "_USE_BOX_" + box + "_monoLatticeFromCubes_Box-" + whichBox + "_numUnits-" + numUnits;
		//prefix = "MD-" + mdCycles + "_firstShellOverwrite_96Options_numUnits-" + numUnits + "_enerCutoff-" + enerCutoff;
		//prefix = "secondShellOverwrite_96Options_" + numUnits + "_";
		//prefix = "MonoBoxesOverwrite_" + whichBox + "_" + numUnits + "_";;
		//prefix = "firstAndSecondShellOverwrite_96Options_" + "%1_" + percentFirst + "_" +  "%2_" + percentSecond + "_" + numUnits + "_";
		//prefix = "firstShellOverwrite_96Options_IncorrectPotential";
		//prefix = "firstShellOverwrite_96Options_LJPotential";
		//prefix = "ZnCl2-3eq_100%RandomOrientation_max3angstromshift";
		//prefix = "ZnCl2-3eq_idealizedCsClBuild_axes=";
		prefix += hydrateAxes.i + "x" + hydrateAxes.j + "x" + hydrateAxes.k + "_";
		prefix += walkLength + "x" + numWalks + "_";
		//prefix = "tetragonalParacrystallineTest";
		//prefix = "FCC_3fold_2fold";
		//prefix = "firstShellShiftedOverwrite_96Options_HyperbolicPotential";
		//prefix = "LJ_firstShellShiftedOverwrite_96Options_";
		//prefix = "LJ_firstShellShiftedOverwrite_110domain_oneDiffractionAxis_";
		//prefix = "MonteCarloTetragonalBuild_distCutoff-" + distCutoff + "_MCcycles-" + monteCarloOrientations;
		//prefix = "LJ_secondShellShiftedOverwrite_96Options_";
		//prefix = "LJ-DEPTH-TEST_firstShellShiftedOverwrite_96Options_" + ljDepth + "-ljDepth_" + numUnits + "-numUnits_";
		//prefix = "LJ-DEPTH-TEST_tetragonal_" + ljDepth + "-ljDepth_" + numUnits + "-numUnits_";
		fName001 = prefix + "_001" + "_" + walkLength + "_" + numWalks + "+" + "(" + numTimes + ").xray";
		fName110 = prefix + "_110" + "_" + walkLength + "_" + numWalks + "+" + "(" + numTimes + ").xray";
		fName111 = prefix + "_111" + "_" + walkLength + "_" + numWalks + "+" + "(" + numTimes + ").xray";
		fName1D = prefix + "_1D" + "_" + walkLength + "_" + numWalks + "+" + "(" + numTimes + ").xray";
		if(numTimes < 0) { numTimes = Integer.MAX_VALUE; }
		ps.println(fName111);
		for(int i = 0; i < numTimes; i++) {
			System.out.println(i);
			c = new Clock();
			c1 = new Clock();
			c.start();
			c1.start();
			mols = l.makeTetragonalFCCLattice();
			//tetra = l.makeTetragonalFCCLattice_withShift();
			//System.out.println("clocks started");
			//tetra = l.alignAlongA110(a);
			//tetra = l.makeTetragonalFCCLattice();
			//tetra = l.makeMonoLatticeFromCubes(a, whichBox);
			//tetra = l.makeMonoLatticeWTetraSpacers(a, cellSize, spacer);
			//tetra = l.monoLatticeRandomLocations(a, whichBox, whichBox, monteCarloOrientations, monteCarloShifts, percentMono);
			//mols = l.makeFirstShellLatticeRandomly(monteCarloOrientations);
			//tetra = l.makeFirstShellLatticeW110Orientation(monteCarloOrientations);
			//tetra = l.makeSecondShellLatticeRandomly(monteCarloOrientations);
			//tetra = l.makeMonoBoxLatticeRandomly(monteCarloOrientations, whichBox, box);
			//tetra = l.makeFirstAndSecondShellLatticeRandomly(monteCarloOrientations, percentFirst, percentSecond);
			//tetra = l.monoBoxesRandomLocations(whichBox, monteCarloOrientations);
			//tetra = l.makeFCC_3fold_2fold();
			//mols = l.createParacrystallineCBr4();
			// mols = l.create3eqZnCl2Lattice(numUnits, numUnits, numUnits, hydrateAxes, false);
			//tetra = l.figure9(monteCarloCycles, idx);
			//tetra = l.monteCarloWalk(monteCarloCycles, idx, numWalks, walkLength);
			//tetra = l.monoLatticeRandomLocationsWRandomSizes(a, cellSize);
			mols = l.monteCarloOrientation(3, monteCarloOrientations);
			//System.out.println("aligned along a 110");
			//tetra = l.dim1Todim3(l.makeTetrahedronLatticeFromXYZ(new File("bigMonoCube.xyz")));
			printLatticeXYZ(molToAtoms(mols, new JVector(1, 1, 1)), i+prefix);
			System.out.println("lattice output before simul");
			//System.out.println("lattice output to file before walk");
			s.setLattice3(mols);
			System.out.println("set lattice");
			//s.calculateTotalEnergy();
			ps.println("ener before: " + s.getTotalEnergy());
			s.walk(walkLength, numWalks, enerCutoff, false, disorderBetweenCycles);
			// there seems to be a problem with the s.MD(walklength) call... Something to do with making all the positions of one molecule null
			//s.MD(walkLength);
			System.out.println("walk complete");
			//s.firstShellWalk(walkLength, 1);
			//s.calculateTotalEnergy();
			ps.println("ener after walk: " + s.getTotalEnergy());
			//walks = s.removeContacts(cutoff, walkLength);
			//s.MD(mdCycles);
			//System.out.println("\"MD\" complete");
			//s.calculateTotalEnergy(); 
			//System.out.println("!\tener after MD: " + s.getTotalEnergy());
			//s.iterate(numWalks, 6);
			c1.stop();
			System.out.println("time to build and simul: " + c1.time() + " ms");
			lattice = molToAtoms(mols, hydrateAxes);
			//lattice = l.create3eqZnCl2();
			//cartesianToFractional(lattice, a);
			//l.rotateLatticeRandomly(lattice, new JVector(numUnits/2*a, numUnits/2*a, numUnits/2*a));
			//lattice = l.getLattice();
			//printLatticeXYZ(molToAtoms(mols, 1), i+prefix);
			printLatticeXYZ(lattice, i+prefix);
			//printLatticeXYZ(lattice, prefix + "outAfter.xyz");
			System.out.println("lattice output to file after walk");
			mols = null;
			//s = null;
			//System.exit(1);
			if(!inEclipse) {
				if(b100) {
					System.out.println("\n\nDiffraction: 001\n\n");
					JavaToC.calcDiffraction(lattice, p001, elemTypes);
					//ewaldProjection(p001);
				}
				if(b110) {
					System.out.println("\n\nDiffraction: 110\n\n");
					JavaToC.calcDiffraction(lattice, p110, elemTypes);
				}
				if(b111) {
					System.out.println("\n\nDiffraction: 111\n\n");
					JavaToC.calcDiffraction(lattice, p111, elemTypes);
				}
				if(bAll){// && i/10 > 0 && i%(i/10) == 0) {
					int idx = 0;
					System.out.println("\n\nDiffraction: 111\n\n");
					// sort the calculated pixels into their correct bins
					JPixel[] tempPix;
					double beginningTime = System.currentTimeMillis();
					double startLoopTime;
					double time;
					for(int j = 0; j < pAll.length; j++) {
						startLoopTime = System.currentTimeMillis();
						System.out.print("idx: " + j + "\tnumPix: " + pAll[j].length);
						tempPix = pAll[j];
						if(tempPix.length > 0) {
							JavaToC.calcDiffraction(lattice, tempPix, elemTypes);
							for(int k = 0; k < tempPix.length; k++) {
								qI[idx++][1] += tempPix[k].getI() * lorentzPolarizationFactor(tempPix[k].getQ().length());
							}
						}
						time = System.currentTimeMillis() - startLoopTime;
						System.out.println("\t" + (time) + "ms for this loop\t" + pAll[j].length/time + " pixels/ms");
					}
					System.out.println("\t" + (System.currentTimeMillis() - beginningTime) + "ms for total calculation of 1d diffraction");
					// coalesce pixels
					
					
				}
			}
			
			addToXYI();
			if(b100) { printPixelsToFile(xyI001, fName001, i+1); }
			if(b110) { printPixelsToFile(xyI110, fName110, i+1); }
			if(b111) { printPixelsToFile(xyI111, fName111, i+1); }
			if(bAll) { printPixelsToFile(coalesce(qI), fName1D, i+1); }
			c.stop();
			System.out.println("\ntime taken for iteration " + i + ": " + c.time() + " ms");
			lattice = null;
			System.gc();
		}
		l = null;
		s = null;
		return walks;
	}
	private double lorentzPolarizationFactor(double Q) {
		double theta = Math.asin(Q*wavelength/4/Math.PI);
		double lorentzPol = (1+Math.pow((Math.cos(2*theta)), 2))/(Math.pow(Math.sin(theta), 2)*Math.cos(theta));
		return lorentzPol;
	}
	public JPixel[] initPixels(JVector[][] vFamily, double qxMax, double qyMax) {
		int numPixels = 0;
		int pixIdx;
		JVector Q, vx, vy;
		JComplex sf;

		for(int i = 0; i < vFamily.length; i++) {
			for(double qy = -qyMax; qy <= qyMax; qy += qStep) {
				for(double qx = -qxMax; qx <= qxMax; qx += qStep) {
					numPixels++;
				}
			}
		}
		JPixel[] p = new JPixel[numPixels];
		System.out.println("numPixels: " + numPixels);
		ComplexScatteringFactor csf = new ComplexScatteringFactor();
		pixIdx = 0;
		for(pixIdx = 0; pixIdx < p.length; pixIdx++) { p[pixIdx] = new JPixel(elemTypes.length); }
		pixIdx = 0;
		/* init new pixels */
		for(int elemTypeIdx = 0; elemTypeIdx < elemTypes.length; elemTypeIdx++) {
			pixIdx = 0;
			for(int i = 0; i < vFamily.length; i++) {
				vx = vFamily[i][1];
				vy = vFamily[i][2];
				sfMap = csf.buildHashMap(Math.max(qxMax, qyMax), qStep, wavelength, elemTypes[elemTypeIdx], vx, vy);
				for(double qy = -qyMax; qy <= qyMax; qy += qStep) {
					for(double qx = -qxMax; qx <= qxMax; qx += qStep) {
						Q = JVector.add(JVector.multiply(vx, qx), JVector.multiply(vy, qy));
						sf = sfMap.get(Q.length());
						
						if(sf == null) {
							try { 
								sf = csf.generateF0(Q, elemTypes[elemTypeIdx], csf.getComplexF(1, wavelength, elemTypes[elemTypeIdx])); 
								//System.out.println("elemTypes[elemTypeIdx]: " + elemTypes[elemTypeIdx]);
							} catch(Exception e) { e.printStackTrace(); }
						}
						//System.out.println(qx + "," + qy+ "," + sf);
						if(elemTypeIdx == 0) { 
							p[pixIdx].setQ((JVector) Q.clone());
						}
						p[pixIdx].setSf(elemTypeIdx, (JComplex) sf.clone());
						p[pixIdx].setI(0);
						pixIdx++;
					}
				}
			}
		}

		return p;
	}
	public JPixel[] initPixels(JVector[] vFamily, double qxMax, double qyMax) {
		int numPixels = 0;
		int pixIdx;
		JVector Q, vx, vy;
		JComplex sf;

		for(double qy = -qyMax; qy <= qyMax; qy += qStep) {
			for(double qx = -qxMax; qx <= qxMax; qx += qStep) {
				numPixels++;
			}
		}
		JPixel[] p = new JPixel[numPixels];
		System.out.println("numPixels: " + numPixels);
		ComplexScatteringFactor csf = new ComplexScatteringFactor();
		pixIdx = 0;
		for(pixIdx = 0; pixIdx < p.length; pixIdx++) { p[pixIdx] = new JPixel(elemTypes.length); }
		pixIdx = 0;
		/* init new pixels */
		for(int elemTypeIdx = 0; elemTypeIdx < elemTypes.length; elemTypeIdx++) {
			pixIdx = 0;
			vx = vFamily[1];
			vy = vFamily[2];
			sfMap = csf.buildHashMap(Math.max(qxMax, qyMax), qStep, wavelength, elemTypes[elemTypeIdx], vx, vy);
			for(double qy = -qyMax; qy <= qyMax; qy += qStep) {
				for(double qx = -qxMax; qx <= qxMax; qx += qStep) {
					Q = JVector.add(JVector.multiply(vx, qx), JVector.multiply(vy, qy));
					sf = sfMap.get(Q.length());
					
					if(sf == null) {
						try { 
							sf = csf.generateF0(Q, elemTypes[elemTypeIdx], csf.getComplexF(1, wavelength, elemTypes[elemTypeIdx])); 
							//System.out.println("elemTypes[elemTypeIdx]: " + elemTypes[elemTypeIdx]);
						} catch(Exception e) { e.printStackTrace(); }
					}
					//System.out.println(qx + "," + qy+ "," + sf);
					if(elemTypeIdx == 0) { 
						p[pixIdx].setQ((JVector) Q.clone());
					}
					p[pixIdx].setSf(elemTypeIdx, (JComplex) sf.clone());
					p[pixIdx].setI(0);
					pixIdx++;
				}
			}
		}

		return p;
	}
	public double[][] initxyI(JVector vx, JVector vy, double qxMax, double qyMax) {
		int numPixels = 0;
		for(double qy = -qyMax; qy <= qyMax; qy += qStep) {
			for(double qx = -qxMax; qx <= qxMax; qx += qStep) {
				numPixels++;
			}
		}
		double[][] xyI = new double[numPixels][3];
		int pixIdx = 0;
		for(double qy = -qyMax; qy <= qyMax; qy += qStep) {
			for(double qx = -qxMax; qx <= qxMax; qx += qStep, pixIdx++) {
				xyI[pixIdx][0] = Math.rint(qx / qStep);
				xyI[pixIdx][1] = Math.rint(qy / qStep);
				xyI[pixIdx][2] = 0;
			}
		}
		System.out.println("xyI done");
		return xyI;
	}
	public double roundQMaxToQStep(double toNorm) {
		return ((int) (toNorm/qStep) ) * qStep;
	}
	public void freeMemory() {
		if(b100) { p001 = null; xyI001 = null; }
		if(b110) { p110 = null; xyI110 = null; }
		if(b111) { p111 = null; xyI111 = null; }
		if(bAll) { pAll = null; qI = null; }
	}
	public void initPixels() {
		JVector vx, vy;
		double qxMax, qyMax;
		qxMax = this.qxMax;
		if(b111) {
			vx = new JVector(1, -1, 0);
			vy = new JVector(-1, -1, 2);
			//qxMax = 1;
			qyMax = qxMax * vx.length()/vy.length();
			//qxMax = 7;
			//qyMax = 4;
			qxMax = roundQMaxToQStep(qxMax);
			qyMax = roundQMaxToQStep(qyMax);
			System.out.println("111");
			System.out.println("vx: " + vx);
			System.out.println("vy: " + vy);
			System.out.println("qxMax: " + qxMax);
			System.out.println("qyMax: " + qyMax);
			p111 = initPixels(JVector.axes111U_ZXY, qxMax, qyMax);
			xyI111 = initxyI(vx, vy, qxMax, qyMax);
		}
		if(b100) {
			vx = new JVector(1, 0, 0);
			vy = new JVector(0, 1, 0);
			//qxMax = 1;
			qyMax = this.qxMax;
			qxMax = roundQMaxToQStep(qxMax);
			qyMax = roundQMaxToQStep(qyMax);
			System.out.println("001");
			System.out.println(vx);
			System.out.println(vy);
			System.out.println(qxMax);
			System.out.println(qyMax);
			p001 = initPixels(JVector.axes100U_ZXY[0], qxMax, qyMax);
			xyI001 = initxyI(vx, vy, qxMax, qyMax);
		}
		if(b110) {
			vx = new JVector(0, 0, 1);
			vy = new JVector(1, -1, 0);
//			qxMax *= .75;
			qyMax = this.qxMax*vx.length() / vy.length();
			qxMax = roundQMaxToQStep(this.qxMax);
			qyMax = roundQMaxToQStep(qyMax);
			System.out.println("110");
			System.out.println(vx);
			System.out.println(vy);
			System.out.println(qxMax);
			System.out.println(qyMax);
			p110 = initPixels(JVector.axes110U_ZXY[0], qxMax, qyMax);
			xyI110 = initxyI(vx, vy, qxMax, qyMax);
		}
		if(bAll) {
			qxMax = 5;
			try {
				Vector<Vector<JPixel>> temp = oneD.createAllPixels(elemTypes, qStep, qxMax, wavelength);
				Vector<JPixel> tempPix;
				pAll = new JPixel[temp.size()][];
				for(int i = 0; i < temp.size(); i++) {
					tempPix = temp.get(i);
					pAll[i] = new JPixel[tempPix.size()];
					pAll[i] = tempPix.toArray(pAll[i]);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
			int idx = 0;
			for(int i = 0; i < pAll.length; i++) {
				idx += pAll[i].length;
			}
			qI = new double[idx][2];
			idx = 0;
			for(int i = 0; i < pAll.length; i++) {
				for(int j = 0; j < pAll[i].length; j++) {
					qI[idx++][0] = pAll[i][j].getQ().length();
				}
			}
		}
	}
	public void readFilesFromFolder(File folder) {
		b100 = false;
		b110 = false;
		b111 = true;
		bAll = false;
		
		wavelength = .13702;
		qStep = 1./10.;
		elemTypes = new int[] {6, 35};
		qxMax = 10;
		double a = 1;
		initPixels();
		File[] filesInFolder = folder.listFiles();
		String outFolder = "D:\\$research\\current\\eclipse projects\\old_JNI_files\\diffraction\\old version of xray calculation";
		outFolder = folder.getParent() + File.separator + "diffraction";
		for(File f : filesInFolder) {
			lattice = XYZTools.readXYZCoordinates(f, a);
			calcDiffractionWSimul(outFolder + File.separator + f.getName());
		}
	}
	/* main method */
	public static void main(String[] args) {
		Diffraction d = new Diffraction();
		File inFolder = new File("D:\\$research\\current\\eclipse projects\\old_JNI_files\\torque_trans-0.5_0.5_Tetragonal_LJPotential_\\xyz");
		inFolder = new File("Z:\\Simulation\\Eric\\CBr4\\Dill-4\\113\\simulOut\\xyz");
		inFolder = new File("C:\\Users\\Eric\\Desktop\\found it\\xyz\\before");
		inFolder = new File("Z:\\Simulation\\Eric\\CBr4\\Dill-4\\135\\simulOut\\xyz");
		d.readFilesFromFolder(inFolder);
		System.exit(1);
		//d.b = BravaisLatticeFactory.getLattice(BravaisLattice.LatticeType.CUBIC_P, new double[] {3.78}, null);
		Clock c = new Clock();
		Date dStart, dFinish;
		int cellSize, spacer, n;
		cellSize = 1;
		spacer = -1;
		n = 5;
		dStart = new Date();
		c.start();
		d.wavelength = .13702;
		//int numUnits = n*(cellSize + spacer + 1)/2;
		int numUnits = 10;
		double torque = .25;
		double trans = .25;
		d.qStep = (int) 1. / ((double)numUnits);
		System.out.println("qStep; " + d.qStep);
		System.out.println("numUnits: " + numUnits);
		//int walkLength = 256;
		//int numWalks = 2;
		int numImages = 25;
		double percentMono = .50;
		//d.ljp = new LennardJonesPotential(8, 17, 2.405426539, .12); // 8 = O, 17 = Cl --> 3eq ZnCl2 H2O
		d.ljp = new LennardJonesPotential(6, 35, 3.37, 0.19);
		d.potLook = new PotentialLookup(.0001, d.ljp);
		//d.elemTypes = new int[] {8, 17, 30};
		d.elemTypes = new int[] {6, 35};

		int numSf = d.elemTypes.length;
		
		int[] monteCarloShift = {0};
		int[] monteCarloOrientations = {0, 2, 4, 6, 8, 10};
		int[] walkLength = {256};
		double percentSecond = .5;
		double percentFirst = .5;
		int[] numWalks = {0};
		//int[] numWalks = {0};
		double ener = d.potLook.lookupPotential(3);
		double[] ljDepth = {.19};
		double[] distCutoff = {3};
		//System.out.println("ljl: " + d.ljl.lookupPotential(3));
		//System.out.println("ljp: " + d.ljp.calcU(3));
		//double[] enerCutoff = {.5};
		//int[] mdCycles = {1, 2, 4, 6, 18, 32, 64};
		int[] mdCycles = {0};
		double[] enerCutoff = {0};
		/*
		int[] monteCarloShift = {100};
		int[] monteCarloOrientations = {100};
		int[] walkLength = {500};
		*/
		int[] boxes = {3};
		d.totalWalks = new int[enerCutoff.length];
		for(int whichBox = 0; whichBox < boxes.length; whichBox++) {
			for(int shift = 0; shift < monteCarloShift.length; shift++) {
				for(int orient = 0; orient < monteCarloOrientations.length; orient++) {
					for(int cycles = 0; cycles < numWalks.length; cycles++) {
						for(int steps = 0; steps < walkLength.length; steps++) {
							for(int cut = 0; cut < enerCutoff.length; cut++) {
								for(int depth = 0; depth < ljDepth.length; depth++) {
									for(int md = 0; md < mdCycles.length; md++) {
										for(int dist = 0; dist < distCutoff.length; dist++) {
											for(int pos = 0; pos < 1; pos++) {
												d.ljp = null;
												d.potLook = null;
												d.ljp = new LennardJonesPotential(17, 8, 3.118145513, ljDepth[depth]);
												d.potLook = new PotentialLookup(.0001, d.ljp);
												d.initPixels();
												d.totalWalks[cut] = d.calcDiffractionWSimul(numImages, numUnits, walkLength[steps], numWalks[cycles], monteCarloOrientations[orient], monteCarloShift[shift],  boxes[whichBox], percentMono,
														percentFirst, percentSecond, enerCutoff[cut], ljDepth[depth], mdCycles[md], distCutoff[dist], pos, .1, .1);
												c.stop();
												dFinish = new Date();
												System.out.println("Started the simulation on:  " + dStart);
												System.out.println("Finished the simulation on: " + dFinish);
												System.out.println("Total simulation time: " + c.time()/1000 + " seconds");
												d.freeMemory();
											}
										}
									}
								}
							}
						}
					}
				}
				
			}
		}
		System.out.println("total walks per cutoff length for " + numImages + " images: ");
		for(int i = 0; i < d.totalWalks.length; i++) {
			System.out.println("cutoff: " + enerCutoff[i] + "\twalks: " + d.totalWalks[i]);
		}
	}
	/* initializing methods */
	static {
		if(!inEclipse) {
			String prop = "java.library.path";
			String curPath = System.getProperty(prop);
			String newPath = ".\\c&cuda";
			System.out.println(prop + " before setting to: " + newPath + ": " + System.getProperty(prop));
			System.setProperty(prop, curPath + ";" + newPath);
			System.out.println(System.getProperty(prop));
			System.out.println(prop + " after setting to: " + newPath + ": " + System.getProperty(prop));
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
			
			JavaToC.cuInit();
		}
		try {
			fos = new FileOutputStream("log.log");
			ps = new PrintStream(fos);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}	
