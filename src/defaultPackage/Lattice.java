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
import java.util.*;
import java.io.*;

public class Lattice {
	final static int totalFails = 50;
	//final static double kb = 1.3806503e-23; // m*m * kg /s/s /K
	final static double kb = 8.617343e-5; // eV K^-1 
	final static double JtoeV = 1./1.602e-19;
	double T = 348;
	double aConstant;
	double distanceCutoff;
	JAtom[] lattice;
	int numUnits;
	int numMols;
	private Simulate s;
	IdealTetrahedron[][] boxes1, boxes2, boxes3, boxes4, boxes5, boxes6, boxes7, firstShells, secondShells, stars;

	private Random r = new Random();
	
	final static double spacing440 = 1.55917;
	final static double shift440 = .5 * spacing440;

	public final static int FIRST_SHELLS = 8;
	public final static int SECOND_SHELLS = 9;
	public final static int STARS = 10;
	
	public Lattice(int numUnits, double aConstant) {
		this.numUnits = numUnits;
		this.aConstant = aConstant;
	}
	
	public static void printOrientations() {
		IdealTetrahedron cur;
		for(int i = 1; i <= 6; i++) {
			cur = getOrientation(i);
//			cur.translate(new JVector(5*i, 0, 0));
			System.out.println(cur.toStringForXYZ());
		}
	}
	
	public static IdealTetrahedron[] getSix4barOrientations() {
		IdealTetrahedron[] tetra = new IdealTetrahedron[6];
		for(int i = 0; i < tetra.length; i++) {
			tetra[i] = getOrientation(i+1);
		}
		
		return tetra;
	}
	/**
	 * Orientations 1 and 2 have a 2-fold aligned with a.
	 * Orientations 3 and 4 have a 2-fold aligned with b.
	 * Orientations 5 and 6 have a 2-fold aligned with c.
	 * @param orientation
	 * @return
	 */
	public static IdealTetrahedron getOrientation(int orientation)
	{
		double bondLen = 1.91;
		double coordinate = Math.sqrt(bondLen*bondLen/3);
		
		JAtom carbonCenter = new JAtom(6, new JVector(0, 0, 0));
		
		JAtom ligand0 = new JAtom(35, new JVector(coordinate, coordinate, coordinate));

		IdealTetrahedron CBr4 = new IdealTetrahedron(carbonCenter, ligand0);
		
		JVector axis = new JVector(1, 0, 0);
		
		JVector origin = new JVector(0, 0, 0);
		
		switch(orientation)
		{
		case 1:
			CBr4.rotate(axis, origin, 45);
			break;
		case 2:
			CBr4.rotate(axis, origin, -45);
			break;
		case 3:
			axis = new JVector(0, 1, 0);
			CBr4.rotate(axis, origin, 45);
			break;
		case 4:
			axis = new JVector(0, 1, 0);
			CBr4.rotate(axis, origin, -45);
			break;
		case 5:
			axis = new JVector(0, 0, 1);
			CBr4.rotate(axis, origin, 45);
			break;
		case 6:
			axis = new JVector(0, 0, 1);
			CBr4.rotate(axis, origin, -45);
			break;
		}
		CBr4.orientation = orientation;
		return CBr4;
	}
	private JVector alignPrimaryBromine_2(IdealTetrahedron[][][] lattice, int a, int b, int c) {
		JVector[] v110s = JVector.v110s;
		
		int[] choices;
		int choice;
		choices = getChoices(lattice, a, b, c);
		choice = choices[r.nextInt(choices.length)];
		lattice[a][b][c].alignBr(v110s[choice]);

		return v110s[choice];
	}
	private void fillInSpaces(IdealTetrahedron[][][] lattice, double aConstant, boolean monteCarloWalk, int monteCarloOrientations, int monteCarloShifts) {
		JVector move;
		
		int numEmpty = 0;
		// get the vector positions of the empty fcc sites
		for(int a = 0; a < numUnits*2; a++) {
			for(int b = 0; b < numUnits*2; b++) {
				for(int c = 0; c < numUnits*2; c++) {
					if((a+b+c) %2 == 0) {
						if(lattice[a][b][c] != null) { continue; }
						numEmpty++;
					}
				}
			}
		}
		JVector[] positions = new JVector[numEmpty];
		numEmpty = 0;
		for(int a = 0; a < numUnits*2; a++) {
			for(int b = 0; b < numUnits*2; b++) {
				for(int c = 0; c < numUnits*2; c++) {
					if((a+b+c) %2 == 0) {
						if(lattice[a][b][c] != null) { continue; }
						positions[numEmpty] = new JVector(a, b, c);
						numEmpty++;
					}
				}
			}
		}
		int a, b, c;
		for(int i = 0; i < positions.length; i++) {
			a = (int) (Math.rint(positions[i].i));
			b = (int) (Math.rint(positions[i].j));
			c = (int) (Math.rint(positions[i].k));
			move = JVector.multiply(new JVector(a, b, c), aConstant/2);
			lattice[a][b][c] = getOrientation((new Random()).nextInt(6));
			lattice[a][b][c].translate(move);
		}
//		monteCarloWalk(lattice, positions, monteCarloOrientations, 0);
//		monteCarloShift(lattice, positions, monteCarloShifts, 0);
	}
	public static Molecule fillInSpaces(Molecule[][][] lattice, double aConstant) {
		JVector move;
		
		int numEmpty = 0;
		// get the vector positions of the empty fcc sites
		for(int a = 0; a < lattice.length; a++) {
			for(int b = 0; b < lattice[a].length; b++) {
				for(int c = 0; c < lattice[a][b].length; c++) {
					if((a+b+c) %2 == 0) {
						if(lattice[a][b][c] != null) { continue; }
						numEmpty++;
					}
				}
			}
		}
		JVector[] positions = new JVector[numEmpty];
		numEmpty = 0;
		for(int a = 0; a < lattice.length; a++) {
			for(int b = 0; b < lattice[a].length; b++) {
				for(int c = 0; c < lattice[a][b].length; c++) {
					if((a+b+c) %2 == 0) {
						if(lattice[a][b][c] != null) { continue; }
						positions[numEmpty] = new JVector(a, b, c);
						numEmpty++;
					}
				}
			}
		}
		int a=0, b=0, c=0;
		for(int i = 0; i < positions.length; i++) {
			a = (int) (Math.rint(positions[i].i));
			b = (int) (Math.rint(positions[i].j));
			c = (int) (Math.rint(positions[i].k));
			move = JVector.multiply(new JVector(a, b, c), aConstant/2);
			lattice[a][b][c] = getOrientation((new Random()).nextInt(6));
			lattice[a][b][c].translate(move);
		}
		return lattice[a][b][c];
	}
	private Molecule randomRotation(Molecule m) {
		
		JVector axis = new JVector(r.nextDouble(), r.nextDouble(), r.nextDouble());
		JVector origin = m.getCenter().getPosition();
		double phi = r.nextDouble()*180;
		m.rotate(axis, origin, phi);
		return m;
	}
	public Molecule[][][] createParacrystallineCBr4() {
		double distortion = .025;
		JVector a1 = new JVector(distortion, distortion, distortion);
		JVector a2 = new JVector(distortion, distortion, distortion);
		JVector a3 = new JVector(distortion, distortion, distortion);
		Vector<Molecule> prevAdditions = new Vector<Molecule>();
		Vector<Molecule> newAdditions = new Vector<Molecule>();
		Molecule[][][] lattice = new Molecule[numUnits*2][numUnits*2][numUnits*2];
		
		lattice[0][0][0] = getOrientation(r.nextInt(6));
		prevAdditions.add(lattice[0][0][0]);
		JVector curPos, newPos, distortionPos;
		Molecule cur;
		int loopIdx = 0;
		while(true) {
			if(prevAdditions.isEmpty()) { break; }
			while(!prevAdditions.isEmpty()) {
				cur = prevAdditions.remove(0);
				curPos = JVector.multiply(cur.getCenter().getPosition(), 2./aConstant);
				curPos.i = Math.rint(curPos.i);
				curPos.j = Math.rint(curPos.j);
				curPos.k = Math.rint(curPos.k);
				for(int i = 0; i < JVector.v110s.length; i++) {
					newPos = JVector.add(JVector.v110s[i], curPos);
					if(withinBounds(newPos, lattice.length) && lattice[(int)newPos.i][(int)newPos.j][(int)newPos.k] == null) {
						newPos = JVector.add(JVector.multiply(JVector.v110s[i], aConstant/2.), cur.getCenter().getPosition());
						//System.out.println(newPos);
						distortionPos = calcParacrystallineDistortion(a1, a2, a3);
						cur = getOrientation(r.nextInt(6));
						cur.translate(JVector.add(newPos, distortionPos));
						newAdditions.add(cur);
						newPos = JVector.multiply(cur.getCenter().getPosition(), 2./aConstant);
						
						System.out.println(newPos);
						newPos.i = Math.rint(newPos.i);
						newPos.j = Math.rint(newPos.j);
						newPos.k = Math.rint(newPos.k);
						if(withinBounds(newPos, lattice.length)) {
							lattice[(int)newPos.i][(int)newPos.j][(int)newPos.k] = cur;
						}
					}
				}
			}
			loopIdx++;
			if(loopIdx % 100 == 0) {
				System.out.println(loopIdx);
			}
			prevAdditions = newAdditions;
			newAdditions = new Vector<Molecule>();
		}
		fillInSpaces(lattice, aConstant); 
		return lattice;
	}
	/**
	 * 
	 * @param a1	distortion parameters along a1
	 * @param a2	distortion parameters along a2
	 * @param a3	distortion parameters along a3
	 * @return
	 */
	private JVector calcParacrystallineDistortion(JVector a1, JVector a2, JVector a3) {
		
		JVector v1 = new JVector(a1.i * r.nextDouble(), a1.j * r.nextDouble(), a1.k * r.nextDouble());
		JVector v2 = new JVector(a2.i * r.nextDouble(), a2.j * r.nextDouble(), a2.k * r.nextDouble());
		JVector v3 = new JVector(a3.i * r.nextDouble(), a3.j * r.nextDouble(), a3.k * r.nextDouble());
		return JVector.add(v1, JVector.add(v2, v3));
	}
	private boolean withinBounds(JVector toTest, int dim) {
		if(toTest.i < 0 || toTest.j < 0 || toTest.k < 0 || toTest.i >= dim || toTest.j >= dim || toTest.k >= dim) { return false; }
		return true;
	}
	public Molecule[][][] create3eqZnCl2Lattice(int a, int b, int c, JVector axes, boolean randomlyRotate) {
		Molecule[] ZnCl4Basis = ZnCl2Hydrate_3eq.getBasisMolecules();
		Molecule[][][] lattice = new Molecule[2*a][2*b][2*c];
		// fill in corners with Zn(OH2)6 molecules
		JVector translate = new JVector();
		
		double maxDisorderAngle = 15;
		for(int i = 0; i < 2*a; i+=2) {
			for(int j = 0; j < 2*b; j+=2) {
				for(int k = 0; k < 2*c; k+=2) {
					lattice[i][j][k] = (IdealOctahedron) ZnCl4Basis[0].clone();
					if(randomlyRotate) {
						lattice[i][j][k].rotate(new JVector(r.nextDouble(), r.nextDouble(), r.nextDouble()), JVector.zero, r.nextDouble()*maxDisorderAngle);
					}
					translate.i = i*axes.i;
					translate.j = j*axes.j;
					translate.k = k*axes.k;
					lattice[i][j][k].translate(translate);
				}
			}
		}
		// fill in centers with ZnCl4 molecules
		boolean orientation = true;
		JVector offset = new JVector();
		for(int i = 1; i < 2*a; i+=2) {
			for(int j = 1; j < 2*b; j+=2) {
				for(int k = 1; k < 2*c; k+=2) {
					if(r.nextInt(2) == 1) {
						lattice[i][j][k] = (IdealTetrahedron) ZnCl4Basis[1].clone();
						orientation = !orientation;
					} else {
						lattice[i][j][k] = (IdealTetrahedron) ZnCl4Basis[2].clone();
						orientation = !orientation;
					}
					if(randomlyRotate) {
						lattice[i][j][k].rotate(new JVector(r.nextDouble(), r.nextDouble(), r.nextDouble()), JVector.zero, r.nextDouble()*maxDisorderAngle);
					}
					translate.i = i*axes.i;
					translate.j = j*axes.j;
					translate.k = k*axes.k;
					lattice[i][j][k].translate(translate);
				}
			}
		}
		return lattice;
	}
	public JAtom[] create3eqZnCl2() {
		JVector[] pos = JVector.getPrimitivePositions(numUnits, 1);
		JVector bccOffset = new JVector(.5, .5, .5);
		JVector newCubeCorner;
		Vector<Molecule> theMolecules = new Vector<Molecule>();
		double znclBondLen = 2.2;
		double znoBondLen = 2.1;
		double clPos = znclBondLen/Math.sqrt(3);
		Molecule zncl4 = new IdealTetrahedron(new JAtom(30, new JVector(0, 0, 0)), new JAtom(17, new JVector(clPos, clPos, clPos)));
		Molecule zno6 = new IdealOctahedron(new JAtom(30, new JVector(0, 0, 0)), 8, znoBondLen);
		Molecule newClPos, newOPos;
		JVector shift = new JVector();
		
		for(int i = 0; i < pos.length; i++) {
			newCubeCorner = JVector.multiply(pos[i], aConstant);
			newClPos = (Molecule) zncl4.clone();
			newClPos.translate(newCubeCorner);
			shift.i = r.nextDouble();
			shift.j = r.nextDouble();
			shift.k = r.nextDouble();
			shift = shift.unit();
			newClPos.translate(shift);
			theMolecules.add(randomRotation(newClPos));
			
			newOPos = (Molecule) zno6.clone();
			newOPos.translate(newCubeCorner);
			newOPos.translate(JVector.multiply(bccOffset, aConstant));

			shift.i = r.nextDouble();
			shift.j = r.nextDouble();
			shift.k = r.nextDouble();
			shift = shift.unit();
			newOPos.translate(shift);
			theMolecules.add(randomRotation(newOPos));
		}
		Vector<JAtom> theAtoms = new Vector<JAtom>();
		Molecule m;
		JAtom[] ligands;
		for(int i = 0; i < theMolecules.size(); i++) {
			m = theMolecules.get(i);
			ligands = m.getLigands();
			theAtoms.add((JAtom) m.getCenter().clone());
			for(int j = 0; j < ligands.length; j++) {
				theAtoms.add((JAtom) ligands[j].clone());
			}
			m = null;
		}
		
		JAtom[] allAtoms = new JAtom[theAtoms.size()];
		allAtoms = theAtoms.toArray(allAtoms);
		return allAtoms;
	}
	public JVector[] getFCCSites() {
		int idx = 0;
		for(int a = 0; a < numUnits*2; a++)
		{
			for(int b = 0; b < numUnits*2; b++)
			{
				for(int c = 0; c < numUnits*2; c++)
				{
					if((a+b+c) %2 == 0) {
						idx++;
					}
				}
			}
		}
		JVector[] sites = new JVector[idx];
		idx = 0;
		for(int a = 0; a < numUnits*2; a++)
		{
			for(int b = 0; b < numUnits*2; b++)
			{
				for(int c = 0; c < numUnits*2; c++)
				{
					if((a+b+c) %2 == 0) {
						sites[idx] = new JVector(a, b, c);
						idx++;
					}
				}
			}
		}
		return sites;
	}
	public boolean canGoHere(JVector pos, int cellSize, IdealTetrahedron[][][] lattice, JVector[] removedSites, JVector[] sites) {
		if(cellSize == 4) 
			System.out.println("boxSize = 4");
		JVector maxLen = JVector.add(pos, new JVector(cellSize*2, cellSize*2, cellSize*2));
		JVector[] fcc = JVector.getFCCPositions(cellSize);
		for(int i = 0; i < fcc.length; i++) {
			fcc[i] = JVector.multiply(fcc[i], 2);
		}
 		if(maxLen.i >= numUnits*2 || maxLen.j >= numUnits*2 || maxLen.k >= numUnits*2) { return false; }

		JVector cur;
		for(int i = 0; i < fcc.length; i++) {
			cur = JVector.add(fcc[i], pos);
			for(int j = 0; j < removedSites.length; j++) {
				if(cur.i == sites[j].i && cur.j == sites[j].j && cur.k == sites[j].k) {
					if(removedSites[j] == null) { return false; }
				}
			}
		}
		return true;
	}
	public IdealTetrahedron[][][] monoLatticeRandomLocationsWRandomSizes(double aConstant, int cellSize) {
		JVector[] sites, removedSites, fcc, tempFCC;
		JVector atom, site, loc;
		boolean insertBox, successfulInsert, random, sequential, tetragonal;
		int failCounter, latticeSelection, boxChoice, successCounter, molCounter, removedCounter, boxSelection,
			totalMols, maxMols;
		IdealTetrahedron[][][] lattice;
		TwoNodeVector tnv;
		DoubleLinkedListVector pos;
		double dist;
		IdealTetrahedron[][] boxes;
		
		pos = new DoubleLinkedListVector(10);
		lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		sites = getFCCSites();
		totalMols = sites.length;
		maxMols = totalMols/2;
		removedSites = getFCCSites();
		
		insertBox = true;
		successfulInsert = true;
		failCounter = 0;
		successCounter = 0;
		molCounter = 0;
		removedCounter = 0;
		random = true;
		sequential = true;
		tetragonal = true;
		if(random) {
			while(insertBox) {
				boxSelection = r.nextInt(3)+2;
				
				switch(boxSelection) {
				case 2: boxes = boxes2; break;
				case 3: boxes = boxes3; break;
				case 4: boxes = boxes4; break;
				default: boxes = boxes2; break;
				}

				fcc = JVector.getFCCPositions(boxSelection);
				tempFCC = new JVector[fcc.length];
				
				latticeSelection = r.nextInt(sites.length);
				if(canGoHere(sites[latticeSelection], cellSize, lattice, removedSites, sites)) {
					
					removedCounter = 0;
					// init the tempFCC array of fcc positions
					for(int j = 0; j < fcc.length; j++) {
						tempFCC[j] = (JVector) fcc[j].clone();
					}
					
					// loop through the molecules
					boxChoice = r.nextInt(boxes.length);
					for(int j = 0; j < boxes[boxChoice].length; j++) {
						// get the atom center
						atom = boxes[boxChoice][j].getCenter().getPosition();
						// loop through the fcc sites
						for(int i = 0; i < fcc.length; i++) {
							if(tempFCC[i] == null) continue; 
							site = JVector.multiply(fcc[i], 8.5);
							// calc the distance from the exact fcc site to the current position
							dist = Math.abs(JVector.subtract(atom, site).length());
							pos.insert(new TwoNodeVector(dist, i, null, null));
						}
						// remove the closest
						tnv = pos.removeHead();
						if(tnv.value < 3) {
							// add the fcc site to the current lattice spot
							loc = JVector.add(JVector.multiply(tempFCC[tnv.key], 2), sites[latticeSelection]).roundInt();
					 		if(loc.i >= numUnits*2 || loc.j >= numUnits*2 || loc.k >= numUnits*2) { continue; }
							lattice[(int) loc.i][(int) loc.j][(int) loc.k] = (IdealTetrahedron) boxes[boxChoice][j].clone();
							lattice[(int) loc.i][(int) loc.j][(int) loc.k].translate(JVector.multiply(sites[latticeSelection], aConstant/2));
							molCounter++;
							if(molCounter == maxMols) { 
								insertBox = false;
								j = boxes[boxChoice].length;
								break;
							}
							// find the JVector in the removedSites array and set it to null
							for(int i = 0; i < removedSites.length; i++) {
								if(loc.i == sites[i].i && loc.j == sites[i].j && loc.k == sites[i].k) {
									removedSites[i] = null;	
									removedCounter++;
								}
							}
							tempFCC[tnv.key] = null;
						}
						pos.clear();
					}
					successfulInsert = true;
				} else { successfulInsert = false; }
	
				if(successfulInsert) { 
					failCounter = 0; 
					successCounter++;
				}
				else { failCounter++; }
				if(failCounter > totalFails) { insertBox = false; }
			}
		}
		// loop through the remaining insertion sites and see if any boxes can fit
		if(sequential && molCounter < maxMols) {
			for(int n = 2; n <= 4; n++) {
				
				switch(n) {
				case 1: boxes = boxes1; break;
				case 2: boxes = boxes2; break;
				case 3: boxes = boxes3; break;
				case 4: boxes = boxes4; break;
				case 5: boxes = boxes5; break;
				case 6: boxes = boxes6; break;
				case 7: boxes = boxes7; break;
				default: boxes = boxes2; break;
				}

				fcc = JVector.getFCCPositions(n);
				tempFCC = new JVector[fcc.length];
				
				for(int a = 0; a < sites.length; a++) {
					if(canGoHere(sites[a], cellSize, lattice, removedSites, sites)) {
		
						// init the tempFCC array of fcc positions
						for(int j = 0; j < fcc.length; j++) {
							tempFCC[j] = JVector.multiply(fcc[j], 1);
						}
						
						// loop through the molecules
						boxChoice = r.nextInt(boxes.length);
						for(int j = 0; j < boxes[boxChoice].length; j++) {
							// get the atom center
							atom = boxes[boxChoice][j].getCenter().getPosition();
							// loop through the fcc sites
							for(int i = 0; i < fcc.length; i++) {
								if(tempFCC[i] == null) { continue; }
								site = JVector.multiply(fcc[i], 8.5);
								// calc the distance from the exact fcc site to the current position
								dist = Math.abs(JVector.subtract(atom, site).length());
								pos.insert(new TwoNodeVector(dist, i, null, null));
							}
							// remove the closest
							tnv = pos.removeHead();
							if(tnv.value < 3) {
								// add the fcc site to the current lattice spot
								loc = JVector.add(JVector.multiply(tempFCC[tnv.key], 2), sites[a]).roundInt();
								lattice[(int) loc.i][(int) loc.j][(int) loc.k] = (IdealTetrahedron) boxes[boxChoice][j].clone();
								lattice[(int) loc.i][(int) loc.j][(int) loc.k].translate(JVector.multiply(sites[a], aConstant/2));
								molCounter++;
								if(molCounter == maxMols) {
									a = sites.length;
									break;
								}
								// find the JVector in the removedSites array and set it to null
								for(int i = 0; i < sites.length; i++) {
									if(removedSites[i] == null) { continue; }
									if(loc.i == removedSites[i].i && loc.j == removedSites[i].j && loc.k == removedSites[i].k) {
										removedSites[i] = null;	
									}
								}
								tempFCC[tnv.key] = null;
							}
							pos.clear();
						}
						successCounter++;
					}
				}
			}
			
		}
		if(tetragonal) { fillInSpaces(lattice, aConstant, false, 0, 0); }
		return lattice;
	}
	
	public IdealTetrahedron[][][] monoLatticeRandomLocations(double aConstant, int cellSize, int whichBox, int monteCarloOrientations, int monteCarloShifts, double percentMono) {
		JVector[] sites, removedSites, fcc, tempFCC;
		JVector atom, site, loc;
		boolean insertBox, successfulInsert, random, sequential, tetragonal;
		int failCounter, latticeSelection, boxChoice, successCounter, molCounter, removedCounter;
		IdealTetrahedron[][][] lattice;
		TwoNodeVector tnv;
		DoubleLinkedListVector pos;
		double dist;
		IdealTetrahedron[][] boxes;
		int numBoxes = 0;
		switch(cellSize) {
		case 1: boxes = boxes1; break;
		case 2: boxes = boxes2; break;
		case 3: boxes = boxes3; break;
		case 4: boxes = boxes4; break;
		case 5: boxes = boxes5; break;
		case 6: boxes = boxes6; break;
		case 7: boxes = boxes7; break;
		default: boxes = boxes2; break;
		}
		boxes = boxes2;
		pos = new DoubleLinkedListVector(10);
		fcc = JVector.getFCCPositions(cellSize);
		tempFCC = new JVector[fcc.length];
		lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		double totalSites = Math.pow(lattice.length/2, 3)*4;
		sites = getFCCSites();
		removedSites = getFCCSites();
		
		insertBox = true;
		successfulInsert = true;
		failCounter = 0;
		successCounter = 0;
		molCounter = 0;
		removedCounter = 0;
		random = true;
		sequential = true;
		tetragonal = true;
		if(random) {
			while(insertBox) {
				if(numBoxes * boxes[0].length < percentMono * totalSites) {
					latticeSelection = r.nextInt(sites.length);
					if(canGoHere(sites[latticeSelection], cellSize, lattice, removedSites, sites)) {
						numBoxes++;
						removedCounter = 0;
						// init the tempFCC array of fcc positions
						for(int j = 0; j < fcc.length; j++) {
							tempFCC[j] = (JVector) fcc[j].clone();
						}
						
						// loop through the molecules
						//boxChoice = r.nextInt(boxes.length);
						if(whichBox < 0) { boxChoice = r.nextInt(boxes.length); }
						else { boxChoice = r.nextInt(12) + 12*whichBox; }
						boxChoice = r.nextInt(boxes.length);
						for(int j = 0; j < boxes[boxChoice].length; j++) {
							// get the atom center
							atom = boxes[boxChoice][j].getCenter().getPosition();
							// loop through the fcc sites
							for(int i = 0; i < fcc.length; i++) {
								if(tempFCC[i] == null) { continue; }
								site = JVector.multiply(fcc[i], 8.5);
								// calc the distance from the exact fcc site to the current position
								dist = Math.abs(JVector.subtract(atom, site).length());
								pos.insert(new TwoNodeVector(dist, i, null, null));
							}
							// remove the closest
							tnv = pos.removeHead();
							if(tnv.value < 3) {
								// add the fcc site to the current lattice spot
								loc = JVector.add(JVector.multiply(tempFCC[tnv.key], 2), sites[latticeSelection]).roundInt();
								lattice[(int) loc.i][(int) loc.j][(int) loc.k] = (IdealTetrahedron) boxes[boxChoice][j].clone();
								lattice[(int) loc.i][(int) loc.j][(int) loc.k].translate(JVector.multiply(sites[latticeSelection], aConstant/2));
								molCounter++;
								// find the JVector in the removedSites array and set it to null
								for(int i = 0; i < removedSites.length; i++) {
									if(loc.i == sites[i].i && loc.j == sites[i].j && loc.k == sites[i].k) {
										removedSites[i] = null;	
										removedCounter++;
									}
								}
								tempFCC[tnv.key] = null;
							}
							pos.clear();
						}
						successfulInsert = true;
					} else { successfulInsert = false; }
		
					if(successfulInsert) { 
						failCounter = 0; 
						successCounter++;
					}
					else { failCounter++; }
					if(failCounter > totalFails) { insertBox = false; }
				}
				else { insertBox = false; }
			}
		}
		System.out.println("number of boxes inserted at random locations: " + numBoxes);
		// loop through the remaining insertion sites and see if any boxes can fit
		if(sequential) {
			for(int a = 0; a < sites.length; a++) {
				if(numBoxes * boxes[0].length < percentMono * totalSites) {
					if(canGoHere(sites[a], cellSize, lattice, removedSites, sites)) {
						numBoxes++;
						// init the tempFCC array of fcc positions
						for(int j = 0; j < fcc.length; j++) {
							tempFCC[j] = JVector.multiply(fcc[j], 1);
						}
						
						// loop through the molecules
						boxChoice = r.nextInt(boxes.length);
						for(int j = 0; j < boxes[boxChoice].length; j++) {
							// get the atom center
							atom = boxes[boxChoice][j].getCenter().getPosition();
							// loop through the fcc sites
							for(int i = 0; i < fcc.length; i++) {
								if(tempFCC[i] == null) { continue; }
								site = JVector.multiply(fcc[i], 8.5);
								// calc the distance from the exact fcc site to the current position
								dist = Math.abs(JVector.subtract(atom, site).length());
								pos.insert(new TwoNodeVector(dist, i, null, null));
							}
							// remove the closest
							tnv = pos.removeHead();
							if(tnv.value < 3) {
								// add the fcc site to the current lattice spot
								loc = JVector.add(JVector.multiply(tempFCC[tnv.key], 2), sites[a]).roundInt();
								lattice[(int) loc.i][(int) loc.j][(int) loc.k] = (IdealTetrahedron) boxes[boxChoice][j].clone();
								lattice[(int) loc.i][(int) loc.j][(int) loc.k].translate(JVector.multiply(sites[a], aConstant/2));
								// find the JVector in the removedSites array and set it to null
								for(int i = 0; i < sites.length; i++) {
									if(removedSites[i] == null) { continue; }
									if(loc.i == removedSites[i].i && loc.j == removedSites[i].j && loc.k == removedSites[i].k) {
										removedSites[i] = null;	
									}
								}
								tempFCC[tnv.key] = null;
							}
							pos.clear();
						}
						successCounter++;
					}
				}
				else break;
			}
		}
		System.out.println("number of total boxes inserted: " + numBoxes);
		System.out.println("Percentage of monoclinic boxes: " + numBoxes * boxes[0].length *100 / (Math.pow(lattice.length/2, 3)*4));
		boolean monteCarloWalk = true;
		if(tetragonal) { fillInSpaces(lattice, aConstant, monteCarloWalk, monteCarloOrientations, monteCarloShifts); }
		return lattice;
	}
	public IdealTetrahedron[][] sort(IdealTetrahedron[][] lattice, int whichBox) {
		DoubleLinkedListMolecule dll = new DoubleLinkedListMolecule(lattice[0].length);
		IdealTetrahedron[][] sorted = new IdealTetrahedron[lattice.length][lattice[0].length];
		//sorted = new IdealTetrahedron[lattice.length][lattice[0].length];
		JVector[] sites = new JVector[JVector.fcc3.length];
		for(int i = 0; i < sites.length; i++) {
			sites[i] = JVector.multiply(JVector.fcc3[i], 8.82/2);
		}
		double dist;
		for(int i = 0; i < lattice.length; i++) {
			for(int k = 0; k < sites.length; k++) {
				dll.clear();
				for(int j = 0; j < lattice[0].length; j++) {
					dist = JVector.subtract(lattice[i][j].center.getPosition(), sites[k]).length();
					dist = Math.abs(dist);
					dll.insert(new TwoNodeMolecule(lattice[i][j], dist, null, null));
				}
				sorted[i][k] = (IdealTetrahedron) dll.removeHead().value.clone();
			}
		}
		return sorted;
	}
	public IdealTetrahedron[][][] monoBoxesRandomLocations(int cellSize, int cycles) {
		JVector[] sites, fcc;
		JVector atom, site, loc;
		int boxChoice, molCounter, posChoice;
		IdealTetrahedron[][][] lattice;
		TwoNodeVector tnv;
		DoubleLinkedListVector pos;
		double dist;
		IdealTetrahedron[][] boxes;
		int numBoxes = 0;
		switch(cellSize) {
		case 1: boxes = boxes1; break;
		case 2: boxes = boxes2; break;
		case 3: boxes = boxes3; break;
		case 4: boxes = boxes4; break;
		case 5: boxes = boxes5; break;
		case 6: boxes = boxes6; break;
		case 7: boxes = boxes7; break;
		default: boxes = boxes2; break;
		}
		pos = new DoubleLinkedListVector(10);
		
		fcc = JVector.fcc3;
		lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		sites = getFCCSites();
		
		molCounter = 0;
		for(int a = 0; a < cycles; a++) {
			for(int b = 0; b < sites.length; b++) {
				// init the tempFCC array of fcc positions

				boxChoice = r.nextInt(boxes.length);
				posChoice = r.nextInt(sites.length);
				
				// loop through the molecules
				for(int j = 0; j < boxes[boxChoice].length; j++) {
					// get the atom center
					atom = boxes[boxChoice][j].getCenter().getPosition();
					// loop through the fcc sites
					for(int i = 0; i < fcc.length; i++) {
						site = JVector.multiply(fcc[i], aConstant/2);
						// calc the distance from the exact fcc site to the current position
						dist = Math.abs(JVector.subtract(atom, site).length());
						pos.insert(new TwoNodeVector(dist, i, null, null));
					}
					// remove the closest
					tnv = pos.removeHead();
					// add the fcc site to the current lattice spot
					loc = JVector.add(JVector.multiply(fcc[tnv.key], 2), sites[posChoice]).roundInt();
					loc = checkPos(loc);
					lattice[(int) loc.i][(int) loc.j][(int) loc.k] = (IdealTetrahedron) boxes[boxChoice][j].clone();
					lattice[(int) loc.i][(int) loc.j][(int) loc.k].translate(JVector.subtract(JVector.multiply(loc, aConstant/2), JVector.multiply(fcc[tnv.key], 2)));
					molCounter++;
					pos.clear();
				}
				numBoxes++;
			}
		}
		System.out.println("number of boxes inserted at random locations: " + numBoxes);
		// loop through the remaining insertion sites and see if any boxes can fit
		return lattice;
	}
	public IdealTetrahedron[][][] makeMonoLatticeWTetraSpacers(double aConstant, int cellSize, int spacer) {
		JVector[] fcc = JVector.fcc2x2x2Pos;
		JVector[] tempFCC = new JVector[fcc.length];
		JVector loc;
		DoubleLinkedListVector pos = new DoubleLinkedListVector(10);
		TwoNodeVector tnv;
		IdealTetrahedron[][][] lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		getS().setLattice3(lattice);
		double dist;
		int choice;
		
		JVector atom, site;
		int primitiveSite = spacer + cellSize + 1;
		IdealTetrahedron[][] boxes;
		switch(cellSize) {
		case 1: boxes = boxes1; break;
		case 2: boxes = boxes2; break;
		case 3: boxes = boxes3; break;
		case 4: boxes = boxes4; break;
		case 5: boxes = boxes5; break;
		case 6: boxes = boxes6; break;
		case 7: boxes = boxes7; break;
		default: boxes = boxes2; break;
		}
		
		for(int a = 0; a < numUnits*2; a++)
		{
			for(int b = 0; b < numUnits*2; b++)
			{
				for(int c = 0; c < numUnits*2; c++)
				{
					if( (a % primitiveSite) == 0 && (b % primitiveSite) == 0 && (c % primitiveSite) == 0)
					{	
						JVector move = new JVector(a, b, c);
						move = JVector.multiply(move, aConstant/2);
						choice = r.nextInt(boxes.length);
						// init the tempFCC array of fcc positions
						for(int j = 0; j < fcc.length; j++) {
							tempFCC[j] = (JVector) fcc[j].clone();
						}
						// loop through the molecules
						for(int j = 0; j < boxes[choice].length; j++) {
							// get the atom center
							atom = boxes[choice][j].getCenter().getPosition();
							// loop through the fcc sites
							for(int i = 0; i < fcc.length; i++) {
								if(tempFCC[i] == null) { continue; }
								site = JVector.multiply(fcc[i], 8.5);
								// calc the distance from the exact fcc site to the current position
								dist = Math.abs(JVector.subtract(atom, site).length());
								pos.insert(new TwoNodeVector(dist, i, null, null));
							}
							// remove the closest
							tnv = pos.removeHead();
							if(tnv.value < 3) {
								// add the fcc site to the current lattice spot
								loc = JVector.add(JVector.multiply(tempFCC[tnv.key], 2), new JVector(a, b, c)).roundInt();
								lattice[(int) loc.k][(int) loc.j][(int) loc.i] = (IdealTetrahedron) boxes[choice][j].clone();
								lattice[(int) loc.k][(int) loc.j][(int) loc.i].translate(move);
								tempFCC[tnv.key] = null;
							}
							pos.clear();
						}
					}
				}
			}
		}
		//System.out.println("complete");
		fillInSpaces(lattice, aConstant, false, 0, 0);
		return lattice;
	}
	public JVector checkPos(JVector curPos) {
		if(curPos.i < 0) { curPos.i += numUnits*2; }
		else if(curPos.i >= numUnits * 2) { curPos.i %= (numUnits * 2); }

		if(curPos.j < 0) { curPos.j += numUnits*2; }
		else if(curPos.j >= numUnits * 2) { curPos.j %= (numUnits * 2); }

		if(curPos.k < 0) { curPos.k += numUnits*2; }
		else if(curPos.k >= numUnits * 2) { curPos.k %= (numUnits * 2); }
		
		return curPos;
	}
	public IdealTetrahedron[][][] makeFirstShellLatticeRandomly(int latticeCycles) {
		IdealTetrahedron[][][] lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		JVector[] sites = getFCCSites();
		int siteSelection, shellSelection, a, b, c;
		
		JVector transVector, curPos, curPos2;
		for(int i = 0; i < latticeCycles; i++) {
			for(int j = 0; j < sites.length; j++) {
				siteSelection = r.nextInt(sites.length);
				shellSelection = r.nextInt(firstShells.length);
				for(int k = 0; k < firstShells[shellSelection].length; k++) {
					transVector = JVector.multiply(sites[siteSelection], aConstant/2);
					curPos = JVector.add(JVector.firstShell[k], sites[siteSelection]);
					curPos2 = checkPos(curPos);
					a = (int) (curPos2.i);
					b = (int) (curPos2.j);
					c = (int) (curPos2.k);
					lattice[a][b][c] = (IdealTetrahedron) firstShells[shellSelection][k].clone();
					lattice[a][b][c].translate(transVector);
				}
			}
		}
//		fillInSpaces(lattice, aConstant, false, 0, 0);
		return lattice;
	}
	public IdealTetrahedron[][][] makeFirstShellLatticeW110Orientation(int latticeCycles) {
		IdealTetrahedron[][][] lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		JVector[] sites = getFCCSites();
		int siteSelection, shellSelection, a, b, c;
		int selection110 = 0;
		
		JVector transVector, curPos, curPos2;
		selection110 = r.nextInt(12);
		for(int i = 0; i < latticeCycles; i++) {
			for(int j = 0; j < sites.length; j++) {
				siteSelection = r.nextInt(sites.length);
				shellSelection = r.nextInt(firstShells.length/firstShells[0].length)*firstShells[0].length+selection110;
				//System.out.print("");
				for(int k = 0; k < firstShells[shellSelection].length; k++) {
					transVector = JVector.multiply(sites[siteSelection], aConstant/2);
					curPos = JVector.add(JVector.firstShell[k], sites[siteSelection]);
					curPos2 = checkPos(curPos);
					a = (int) (curPos2.i);
					b = (int) (curPos2.j);
					c = (int) (curPos2.k);
					lattice[a][b][c] = (IdealTetrahedron) firstShells[shellSelection][k].clone();
					lattice[a][b][c].translate(transVector);
				}
			}
		}
		fillInSpaces(lattice, aConstant, false, 0, 0);
		return lattice;
	}
	public IdealTetrahedron[][][] makeSecondShellLatticeRandomly(int latticeCycles) {
		IdealTetrahedron[][][] lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		JVector[] sites = getFCCSites();
		int siteSelection, shellSelection, a, b, c;
		
		JVector transVector, curPos, curPos2;
		for(int i = 0; i < latticeCycles; i++) {
			for(int j = 0; j < sites.length; j++) {
				siteSelection = r.nextInt(sites.length);
				shellSelection = r.nextInt(secondShells.length);
				for(int k = 0; k < secondShells[shellSelection].length; k++) {
					transVector = JVector.multiply(sites[siteSelection], aConstant/2);
					curPos = JVector.add(JVector.secondShell[k], sites[siteSelection]);
					curPos2 = checkPos(curPos);
					a = (int) (curPos2.i);
					b = (int) (curPos2.j);
					c = (int) (curPos2.k);
					lattice[a][b][c] = (IdealTetrahedron) secondShells[shellSelection][k].clone();
					lattice[a][b][c].translate(transVector);
				}
			}
		}
//		fillInSpaces(lattice, aConstant, false, 0, 0);
		return lattice;
	}

	public JVector[] getStarLocations() {
		JVector[] starLocations = new JVector[stars[0].length];
		
		starLocations[0] = new JVector(0, 0, 0);
		
		int scalar;
		int numPerBox = stars[0].length;
		for(int i = 1; i < numPerBox; i++) {
			scalar = i / JVector.firstShell.length + 1;
			starLocations[i] = (JVector) JVector.firstShell[(i-1) % (JVector.firstShell.length - 1) + 1].clone();
			starLocations[i].multiply(scalar);
		}
		
		return starLocations;
	}
	
	public IdealTetrahedron[][][] makeStarLatticeRandomly(int latticeCycles) {
		IdealTetrahedron[][][] lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		JVector[] sites = getFCCSites();
		int siteSelection, starSelection, a, b, c;
		
		JVector[] starLocations = getStarLocations();
		
		JVector transVector, curPos, curPos2;
		for(int i = 0; i < latticeCycles; i++) {
			for(int j = 0; j < sites.length; j++) {
				siteSelection = r.nextInt(sites.length);
				starSelection = r.nextInt(stars.length);
				for(int k = 0; k < stars[starSelection].length; k++) {
					IdealTetrahedron toPlace = stars[starSelection][k];
					if(toPlace == null)
						continue;
					
					transVector = JVector.multiply(sites[siteSelection], aConstant/2);
					curPos = JVector.add(starLocations[k], sites[siteSelection]);
					curPos2 = checkPos(curPos);
					a = (int) (curPos2.i);
					b = (int) (curPos2.j);
					c = (int) (curPos2.k);
					IdealTetrahedron tetra = (IdealTetrahedron) toPlace.clone();
					tetra.translate(transVector);
					JVector difference = JVector.subtract(tetra.getCenter().getPosition(), 
							JVector.multiply(curPos, aConstant/2));
					lattice[a][b][c] = tetra;
				}
			}
		}
		fillInSpaces(lattice, aConstant, false, 0, 0);
		return lattice;
	}
	public IdealTetrahedron[][][] makeFirstAndSecondShellLatticeRandomly(int latticeCycles, double percentFirst, double percentSecond) {
		IdealTetrahedron[][][] lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		JVector[] sites = getFCCSites();
		int siteSelection, shellSelection, a, b, c;
		
		JVector transVector, curPos, curPos2;
		int numFirst = 0, numSecond = 0;
		// scale the percentages
		percentFirst *= ((double) secondShells[0].length) /((double) firstShells[0].length);
		double totalPercent = percentFirst + percentSecond;
		percentFirst /= totalPercent;
		percentSecond /= totalPercent;
		System.out.println("percentFirst: " + percentFirst + "\tpercentSecond: " + percentSecond + "\ttotalPercent: " + (percentFirst + percentSecond));
		double whichShell;
		IdealTetrahedron[][] shell;
		JVector[] shellSites;
		for(int i = 0; i < latticeCycles; i++) {
			for(int j = 0; j < sites.length; j++) {
				siteSelection = r.nextInt(sites.length);
				whichShell = r.nextDouble();
				if(whichShell < percentSecond) {shell = secondShells; shellSites = JVector.secondShell; numSecond++;}
				else { shell = firstShells; shellSites = JVector.firstShell; numFirst++; }
				shellSelection = r.nextInt(shell.length);
				for(int k = 0; k < shell[shellSelection].length; k++) {
					transVector = JVector.multiply(sites[siteSelection], aConstant/2);
					curPos = JVector.add(shellSites[k], sites[siteSelection]);
					curPos2 = checkPos(curPos);
					a = (int) (curPos2.i);
					b = (int) (curPos2.j);
					c = (int) (curPos2.k);
					lattice[a][b][c] = (IdealTetrahedron) shell[shellSelection][k].clone();
					lattice[a][b][c].translate(transVector);
				}
			}
		}
		System.out.println("numFirst: " + numFirst + "\tnumSecond: " + numSecond + "\tnumTotal: " + (numFirst + numSecond));
		System.out.println("percentFirst: " + ((double) numFirst)/((double) numFirst + numSecond) + "\tpercentSecond: " + ((double) numSecond)/((double) numFirst + numSecond));
		fillInSpaces(lattice, aConstant, false, 0, 0);
		return lattice;
	}
	public IdealTetrahedron[][][] makeMonoBoxLatticeRandomly(int latticeCycles, int layers, int box) {
		IdealTetrahedron[][][] lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		JVector[] sites = getFCCSites();
		JVector[] fcc;
		IdealTetrahedron[][] boxes;
		switch(layers) {
		case 1: fcc = JVector.fcc1; boxes = boxes1; break;
		case 2: fcc = JVector.fcc2; boxes = boxes2; break;
		case 3: fcc = JVector.fcc3; boxes = boxes3; break;
		case 4: fcc = JVector.fcc4; boxes = boxes4; break;
		case 5: fcc = JVector.fcc5; boxes = boxes5; break;
		default: fcc = JVector.fcc3; boxes = boxes3; break;
		}
		boxes = sort(boxes, layers);
		int siteSelection, boxSelection, a, b, c;
		
		JVector transVector, curPos, curPos2;
		for(int i = 0; i < latticeCycles; i++) {
			for(int j = 0; j < sites.length; j++) {
				siteSelection = r.nextInt(sites.length);
				boxSelection = r.nextInt(12)+12*box;
				for(int k = 0; k < boxes[boxSelection].length; k++) {
					transVector = JVector.multiply(sites[siteSelection], aConstant/2);
					curPos = JVector.add(fcc[k], sites[siteSelection]);
					curPos2 = checkPos(curPos);
					a = (int) (curPos2.i);
					b = (int) (curPos2.j);
					c = (int) (curPos2.k);
					lattice[a][b][c] = (IdealTetrahedron) boxes[boxSelection][k].clone();
					lattice[a][b][c].translate(transVector);
				}
			}
		}
		fillInSpaces(lattice, aConstant, false, 0, 0);
		return lattice;
	}
	public IdealTetrahedron[][][] makeMonoBoxesLatticeRandomly(int latticeCycles) {
		IdealTetrahedron[][][] lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		JVector[] sites = getFCCSites();
		int siteSelection, shellSelection, a, b, c;
		
		JVector transVector, curPos, curPos2;
		JVector[] fcc = JVector.fcc3;
		boxes2 = sort(boxes2, 2);
		for(int i = 0; i < latticeCycles; i++) {
			for(int j = 0; j < sites.length; j++) {
				siteSelection = r.nextInt(sites.length);
				shellSelection = r.nextInt(boxes2.length);
				for(int k = 0; k < boxes2[shellSelection].length; k++) {
					transVector = JVector.multiply(sites[siteSelection], aConstant/2);
					curPos = JVector.add(fcc[k], sites[siteSelection]);
					curPos2 = checkPos(curPos);
					a = (int) (curPos2.i);
					b = (int) (curPos2.j);
					c = (int) (curPos2.k);
					lattice[a][b][c] = (IdealTetrahedron) boxes2[shellSelection][k].clone();
					lattice[a][b][c].translate(transVector);
				}
			}
		}
		fillInSpaces(lattice, aConstant, false, 0, 0);
		return lattice;
	}
	public IdealTetrahedron[][][] makeMonoLatticeFromCubes(double aConstant, int whichBox)
	{
		
		JVector[] fcc = JVector.fcc2x2x2Pos;
		JVector[] tempFCC = new JVector[fcc.length];
		JVector loc;
		DoubleLinkedListVector pos = new DoubleLinkedListVector(10);
		TwoNodeVector tnv;
		IdealTetrahedron[][][] lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		double dist;
		int choice;
		
		JVector atom, site;
		int maxChoice = boxes3.length/12;
		choice = whichBox;
		JVector[] positions = JVector.getPrimitivePositions(numUnits*2, whichBox+1);
		for(int a = 0; a < positions.length; a++) {
			choice = r.nextInt(maxChoice) + r.nextInt(12);
			//else { choice = r.nextInt(12) + 12*whichBox; }
			//choice = r.nextInt(maxChoice);
			JVector move = JVector.multiply(positions[a], aConstant/2);
			// init the tempFCC array of fcc positions
			for(int j = 0; j < fcc.length; j++) {
				tempFCC[j] = (JVector) fcc[j].clone();
			}
			// loop through the molecules
			for(int j = 0; j < boxes3[choice].length; j++) {
				// get the atom center
				atom = boxes3[choice][j].getCenter().getPosition();
				// loop through the fcc sites
				for(int i = 0; i < fcc.length; i++) {
					if(tempFCC[i] == null) { continue; }
					site = JVector.multiply(fcc[i], 8.5);
					// calc the distance from the exact fcc site to the current position
					dist = Math.abs(JVector.subtract(atom, site).length());
					pos.insert(new TwoNodeVector(dist, i, null, null));
				}
				// remove the closest
				tnv = pos.removeHead();
				if(tnv.value < 3) {
					// add the fcc site to the current lattice spot
					loc = JVector.add(JVector.multiply(tempFCC[tnv.key], 2), positions[a]);
					lattice[(int) loc.k][(int) loc.j][(int) loc.i] = (IdealTetrahedron) boxes3[choice][j].clone();
					lattice[(int) loc.k][(int) loc.j][(int) loc.i].translate(move);
					tempFCC[tnv.key] = null;
				}
				pos.clear();
			}
		}
		//System.out.println("complete");
		return lattice;
	}
	private boolean isAligned(JVector alignment, IdealTetrahedron cur) {
		JVector[] atoms = cur.getVectors();
		double phi = 0;
		for(int i = 1; i < atoms.length; i++) {
			phi = JVector.angle(alignment, atoms[i]);
			if(phi == 0)
				return true;
		}
		return false;
	}
	private int checkCoord(int coord, int boxSize) { return (coord + boxSize) % boxSize; }
	public void rotateLatticeRandomly(JAtom[] lattice, JVector origin) {
		JVector axis = new JVector(Math.random(), Math.random(), Math.random());
		double phi = Math.random()*180;
		
		for(int i = 0; i < lattice.length; i++) {
			lattice[i].setPosition(JVector.rotate(lattice[i].getPosition(), axis, origin, phi));
		}
	}
	private int[] getChoices(IdealTetrahedron[][][] lattice, int a, int b, int c) {
		Stack<Integer> choices = new Stack<Integer>();
		int testA, testB, testC;
		for(int i = 0; i < JVector.v110s.length; i++) {
			testA = a + (int) JVector.v110s[i].i;
			testB = b + (int) JVector.v110s[i].j;
			testC = c + (int) JVector.v110s[i].k;
			testA = checkCoord(testA, lattice.length);
			testB = checkCoord(testB, lattice.length);
			testC = checkCoord(testC, lattice.length);
			if(lattice[testA][testB][testC] == null || !isAligned(JVector.multiply(JVector.v110s[i], -1), lattice[testA][testB][testC])) {
				choices.push(i);
			}
		}
		int[] options = new int[choices.size()];
		for(int i = 0; i < options.length; i++) {
			options[i] = choices.pop();
		}
		return options;
	}
	private void align2Fold(IdealTetrahedron[][][] lattice, int a, int b, int c) {
		JVector[] v110s = JVector.v110s;
		
		int[] choices;
		int choice;
		choices = getChoices(lattice, a, b, c);
		choice = choices[r.nextInt(choices.length)];
		lattice[a][b][c].align2Fold(v110s[choice]);
		try {
			lattice[a][b][c].translate(JVector.multiply(v110s[choice].unit(), .2));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	private JVector alignPrimaryBromine(IdealTetrahedron[][][] lattice, int a, int b, int c) {
		JVector[] v110s = JVector.v110s;
		
		int[] choices;
		int choice;
		choices = getChoices(lattice, a, b, c);
		choice = choices[r.nextInt(choices.length)];
		lattice[a][b][c].alignBr(v110s[choice]);
		try {
			lattice[a][b][c].translate(JVector.multiply(v110s[choice].unit(), .2));
		} catch (Exception e) {
			e.printStackTrace();
		}
		return v110s[choice];
	}
	private void alignSecondaryBromine(IdealTetrahedron[][][] lattice, int a, int b, int c, JVector primaryAlignment) {
		IdealTetrahedron temp;
		DoubleLinkedListVector dllv = new DoubleLinkedListVector(6);
		for(int i = 0; i < 6; i++) {
			temp = null;
			temp = (IdealTetrahedron) lattice[a][b][c].clone();
			temp.align2ndBr(primaryAlignment, i);
			dllv.insert(new TwoNodeVector(getS().calculateEnergyofFirstShell(temp), i, null, null));
		}
		int secondaryAlignment = dllv.removeHead().key;
		lattice[a][b][c].align2ndBr(primaryAlignment, secondaryAlignment);

		try {
			lattice[a][b][c].translate(JVector.multiply(JVector.v110s[secondaryAlignment].unit(), .1));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public IdealTetrahedron[][][] alignAlongA110(double aConstant) {
		IdealTetrahedron[][][] lattice = makeTetragonalFCCLattice();
		getS().setLattice3(lattice);
		JVector[][][] primaryAlignments = new JVector[lattice.length][lattice[0].length][lattice[0][0].length];
		for(int a = 0; a < lattice.length; a++) {
			for(int b = 0; b < lattice[a].length; b++) {
				for(int c = 0; c < lattice[a][b].length; c++) {
					if(lattice[a][b][c] == null) { continue; }
					primaryAlignments[a][b][c] = alignPrimaryBromine(lattice, a, b, c);
				}
			}
		}
		for(int a = 0; a < lattice.length; a++) {
			for(int b = 0; b < lattice[a].length; b++) {
				for(int c = 0; c < lattice[a][b].length; c++) {
					if(lattice[a][b][c] == null) { continue; }
					alignSecondaryBromine(lattice, a, b, c, primaryAlignments[a][b][c]);
				}
			}
		}
		return lattice;
	}
	public IdealTetrahedron[][][] makeRandomlyOrientedFCCLattice() {
		IdealTetrahedron[][][] lattice = makeTetragonalFCCLattice();
		
		double phix, phiy, phiz;
		JVector origin;
		for(IdealTetrahedron[][] tet2 : lattice)
			for(IdealTetrahedron[] tet1 : tet2)
				for(IdealTetrahedron tet : tet1) {
					if(tet == null)
						continue;
					
					phix = r.nextDouble() * 180;
					phiy = r.nextDouble() * 180;
					phiz = r.nextDouble() * 180;
					origin = tet.center.getPosition();
					
					tet.rotate(JVector.add(JVector.x, origin), origin, phix);
					tet.rotate(JVector.add(JVector.y, origin), origin, phiy);
					tet.rotate(JVector.add(JVector.z, origin), origin, phiz);
				}
		
		return lattice;
	}
	public IdealTetrahedron[][][] makeTetragonalFCCLattice()
	{
		IdealTetrahedron[][][] lattice = new IdealTetrahedron[numUnits*2][numUnits*2][numUnits*2];
		for(int c = 0; c < numUnits*2; c++) {
			for(int b = 0; b < numUnits*2; b++) {
				for(int a = 0; a < numUnits*2; a++) {
					if( (a + b + c) % 2 == 0) {						
						IdealTetrahedron temp = getOrientation(r.nextInt(6) + 1);
						//IdealTetrahedron temp = getOrientation((int)( Math.random() * 0) + 1);
						JVector move = new JVector(a, b, c);
						
						move = JVector.multiply(move, aConstant / 2);
						
						temp.translate(move);
						
						lattice[a][b][c] = (IdealTetrahedron) temp.clone();
					}
				}
			}
		}
		//System.out.println("complete");
		return lattice;
	}
	public IdealTetrahedron[][][] makeTetragonalFCCLattice_withShift() {
		IdealTetrahedron[][][] lattice = makeTetragonalFCCLattice();
		getS().setLattice3(lattice);
		JVector trans;
		for(int c = 0; c < numUnits*2; c++) {
			for(int b = 0; b < numUnits*2; b++) {
				for(int a = 0; a < numUnits*2; a++) {
					if( (a + b + c) % 2 == 0) {
						trans = new JVector(0, 0, 0);
						for(int j = 2; j < 4; j++) {
							trans = JVector.add(trans, JVector.subtract(lattice[a][b][c].ligands[j].getPosition(), lattice[a][b][c].center.getPosition()));
						}
						trans = JVector.multiply(trans, 1./3./2.);
						lattice[a][b][c].translate(trans);
					}
				}
			}
		}
		return lattice;
	}
	public IdealTetrahedron[][][] makeFCC_3fold_2fold() {
		
		IdealTetrahedron[][][] lattice = makeTetragonalFCCLattice();
		getS().setLattice3(lattice);
		JVector trans;
		for(int c = 0; c < numUnits*2; c++) {
			for(int b = 0; b < numUnits*2; b++) {
				for(int a = 0; a < numUnits*2; a++) {
					if( (a + b + c) % 2 == 0) {
						trans = new JVector(0, 0, 0);
						switch(r.nextInt(2)) {
						case 0:
							alignSecondaryBromine(lattice, a, b, c, alignPrimaryBromine(lattice, a, b, c));
							for(int j = 1; j < 4; j++) {
								trans = JVector.add(trans, JVector.subtract(lattice[a][b][c].ligands[j].getPosition(), lattice[a][b][c].center.getPosition()));
							}
							trans = JVector.multiply(trans, 1./3./2.);
							lattice[a][b][c].translate(trans);
							break;
						case 1:
							for(int j = 2; j < 4; j++) {
								trans = JVector.add(trans, JVector.subtract(lattice[a][b][c].ligands[j].getPosition(), lattice[a][b][c].center.getPosition()));
							}
							trans = JVector.multiply(trans, 1./3./2.);
							lattice[a][b][c].translate(trans);
							break;
						}
						
					}
				}
			}
		}
		return lattice;
	}
	public void buildRandomLattice(int numScatterers, int scattererType) {
		lattice = new JAtom[numScatterers];
		JVector coord = new JVector();
		for(int i = 0; i < lattice.length; i++) {
			coord.setI(Math.random());
			coord.setJ(Math.random());
			coord.setK(Math.random());
			coord = JVector.multiply(coord, numUnits);
			lattice[i] = new JAtom(scattererType, coord);
		}
	}
	public void build100Lattice(int numScatterers, int scattererType, String xyz) {
		lattice = new JAtom[numScatterers];
		JVector coord = new JVector();
		
		for(int i = 0; i < lattice.length; i++) {
			coord.setI(Math.random());
			coord.setJ(Math.random());
			coord.setK(Math.random());
			coord = JVector.multiply(coord, numUnits);
			if(xyz.contains("x")) {
				coord.setI(r.nextInt(numUnits));
			}
			if(xyz.contains("y")) {
				coord.setJ(r.nextInt(numUnits));
			}
			if(xyz.contains("z")) {
				coord.setK(r.nextInt(numUnits));
			}
			lattice[i] = new JAtom(scattererType, coord);
		}
	}
	public void buildFCCLatticeWDisplace(int scattererType, String xyz, double percent) {
		lattice = new JAtom[(int) (Math.pow(numUnits, 3) * 4)];
		int scattererIdx = 0;
		JVector coord;
		double amt;
		
		for(int a = 0; a < numUnits*2; a++) {
			for(int b = 0; b < numUnits*2; b++) {
				for(int c = 0; c < numUnits*2; c++) {
					if((a + b + c) % 2 == 0) {
						coord = new JVector(a, b, c);
						if(xyz.contains("x")) {
							amt = (((((r.nextDouble() * 2)-1)*percent)));
							coord.setI(amt+a);
						}
						if(xyz.contains("y")) {
							amt = (((((r.nextDouble() * 2)-1)*percent)));
							coord.setJ(amt+b);
						}
						if(xyz.contains("z")) {
							amt = (((((r.nextDouble() * 2)-1)*percent)));
							coord.setK(amt+c);
						}
						if(xyz.contains("xy")) {
							amt = (((((r.nextDouble() * 2)-1)*percent)));
							coord.setI(amt+a);
							coord.setJ(amt+b);
						} else if(xyz.contains("xz")) {
							amt = (((((r.nextDouble() * 2)-1)*percent)));
							coord.setI(amt+a);
							coord.setK(amt+c);
						} else if(xyz.contains("yz")) {
							amt = (((((r.nextDouble() * 2)-1)*percent)));
							coord.setJ(amt+b);
							coord.setK(amt+c);
						}
						coord = JVector.multiply(coord, .5);
						lattice[scattererIdx] = new JAtom(scattererType, coord);
						scattererIdx++;
					}
				}
			}
		}
		System.out.println("num scatterers: " + scattererIdx);
	}
	public IdealTetrahedron[][][] figure7() {
		IdealTetrahedron[][][] lattice = makeTetragonalFCCLattice();
		
		return lattice;
	}
	public IdealTetrahedron[][][] figure9(int cycles, int idx) {
		IdealTetrahedron[][][] lattice = makeTetragonalFCCLattice();
		numMols = (int) (Math.pow(lattice.length/2, 3) * 4);
		getS().setLattice3(lattice);
		double enerBefore, enerAfter;
		int newChoice, oldChoice, twoOrThree;
		int total3fold = 0;
		int threeFoldReset = 0;
		int total2fold = 0;
		int twoFoldReset = 0;
		
		JVector move = new JVector();
		double e;
		for(int i = 0; i < cycles; i++) {
			for(int a = 0; a < lattice.length; a++) {
				for(int b = 0; b < lattice[a].length; b++) {
					for(int c = 0; c < lattice[a][b].length; c++) {
						if(lattice[a][b][c] != null) {
							move.i = lattice[a][b][c].center.getPosition().i;
							move.j = lattice[a][b][c].center.getPosition().j;
							move.k = lattice[a][b][c].center.getPosition().k;
							
							enerBefore = getS().calculateEnergyofFirstShell(lattice[a][b][c]);
							oldChoice = lattice[a][b][c].orientation;
							twoOrThree = r.nextInt(2);
							switch(twoOrThree) {
							// align 2 fold with cubic axes
							case 0:
								total2fold++;
								newChoice = r.nextInt(6);
								while(newChoice == oldChoice) { newChoice = r.nextInt(6); }
								lattice[a][b][c] = getOrientation(newChoice);
								lattice[a][b][c].translate(move);
								break;
							case 1:		
								total3fold++;
								// align 3 fold along intermolecular axis
								JVector alignment = alignPrimaryBromine(lattice, a, b, c);
								if(idx == 1)
									// also align a secondary bromine
									alignSecondaryBromine(lattice, a, b, c, alignment);
								break;
							}
							enerAfter = getS().calculateEnergyofFirstShell(lattice[a][b][c]);
							if(enerAfter >= enerBefore) {
								e = -Math.abs((enerAfter-enerBefore))/ kb/T;
								e = Math.exp(e);
								if(e < r.nextDouble()) {
									lattice[a][b][c] = getOrientation(oldChoice);
									lattice[a][b][c].translate(move);	
									if(twoOrThree == 1) { threeFoldReset++; }
									else { twoFoldReset++; }
										
								}
							}
						}
					}
				}
			}
		}
		System.out.println("Total number of 3-fold alignments: " + total3fold + "\t3-fold reset: " + threeFoldReset);
		System.out.println("Total number of 2-fold alignments: " + total2fold + "\t2-fold reset: " + twoFoldReset);
		return lattice;
	}
	public boolean checkDistances(IdealTetrahedron center, double distanceCutoff) {
		DoubleLinkedListMolecule list = getS().getSurrounding(center);
		Molecule cur;
		int listLen = list.length;
		JAtom[] cent = center.ligands;
		JAtom[] surr;
		for(int i = 0; i < listLen; i++) {
			cur = list.removeHead().value;
			for(int br1 = 0; br1 < cent.length; br1++) {
				surr = cur.getLigands();
				for(int br2 = 0; br2 < surr.length; br2++) {
					if(JVector.subtract(cent[br1].getPosition(), surr[br2].getPosition()).length() < distanceCutoff)
						return true;
				}
			}
		}
		return false;
	}
	public IdealTetrahedron[][][] monteCarloOrientation(double distanceCutoff, int cycles) {
		IdealTetrahedron[][][] lattice = makeTetragonalFCCLattice();
		getS().setLattice3(lattice);
		numMols = (int) (Math.pow(lattice.length/2, 3) * 4);
		getS().setLattice3(lattice);
		double enerBefore, enerAfter;
		int newChoice, oldChoice;
		int total2fold = 0;
		int twoFoldReset = 0;
		int distCutoff = 0;
		
		JVector move = new JVector();
		double e;
		//int walkTime = cycles / latticeCycles;
		boolean newOrientation = false;
		double energy;
		int numTotalUnphysicalContacts = 0;
		int numUnphysicalRemoved = 0;
		getS().calculateTotalEnergy();
		System.out.println("Target lattice energy: " + -1 * numMols * 0.763 + " eV");
		System.out.println("Starting lattice energy: " + getS().getTotalEnergy() + " eV");
		int a, b, c;
		for(int i = 0; i < cycles; i++) {
			/*if(cycles % walkTime == 0) {
				s.walk(walkLength, 1);
			}*/
			Vector<JVector> positions = new Vector<JVector>();
			for(a = 0; a < lattice.length; a++) {
				for(b = 0; b < lattice[a].length; b++) {
					for(c = 0; c < lattice[a][b].length; c++) {
						if(lattice[a][b][c] != null) {
							positions.add(new JVector(a, b, c));
						}
					}
				}
			}
			
			getS().calculateTotalEnergy();
			numTotalUnphysicalContacts = 0;
			numUnphysicalRemoved = 0;
			
			while(positions.size() > 0) {
				JVector curPos = positions.remove(r.nextInt(positions.size()));
				a = (int) curPos.i;
				b = (int) curPos.j;
				c = (int) curPos.k;
				newOrientation = checkDistances(lattice[a][b][c], distanceCutoff);
				total2fold++;
				if(newOrientation) {
					distCutoff++;
					numTotalUnphysicalContacts++;
					enerBefore = getS().calculateEnergyofFirstShell(lattice[a][b][c]);
					oldChoice = lattice[a][b][c].orientation;
					move = lattice[a][b][c].getCenter().getPosition();
					
					do { newChoice = r.nextInt(6); }
					while(newChoice == oldChoice);
					
					lattice[a][b][c] = getOrientation(newChoice);
					
					
					lattice[a][b][c].translate(move);
					
					enerAfter = getS().calculateEnergyofFirstShell(lattice[a][b][c]);
					if(enerAfter > enerBefore) {
						e = -Math.abs((enerAfter-enerBefore))/ kb/T;
						e = Math.exp(e);
						double chance = r.nextDouble();
						if(e > chance) {
							lattice[a][b][c] = getOrientation(oldChoice);
							lattice[a][b][c].translate(move);	
							twoFoldReset++;
								
						}
					} else {
						numUnphysicalRemoved++;
					}
				}
			}
			System.out.print("Total lattice energy after step: " + (i+1) + " of " + cycles + ": " + getS().getTotalEnergy() + " eV");
			System.out.print("\tNumber of unphysical contacts: " + numTotalUnphysicalContacts);
			System.out.println("\tNumber of unphysical contacts removed: " + numUnphysicalRemoved);
		}
		System.out.println("Total number of 2-fold tests: " + total2fold + "\tdistanceCutoff " + distCutoff + "\t2-fold reset: " + twoFoldReset);
		System.out.println("Total number of changed orientations: " + (distCutoff-twoFoldReset));
		return lattice;
	}
	public IdealTetrahedron[][][] monteCarloShift(int cycles) {
		IdealTetrahedron[][][] lattice = makeTetragonalFCCLattice();
		getS().setLattice3(lattice);
		numMols = (int) (Math.pow(lattice.length/2, 3) * 4);
		getS().setLattice3(lattice);
		double enerBefore, enerAfter;
		int total3fold = 0;
		int threeFoldReset = 0;
		int total2fold = 0;
		int twoFoldReset = 0;
		
		JVector translate = new JVector();
		JVector cPos;
		JVector brPos;
		int br1, br2;
		double e;
		//int walkTime = cycles / latticeCycles;
		for(int i = 0; i < cycles; i++) {
			/*if(cycles % walkTime == 0) {
				s.walk(walkLength, 1);
			}*/
			for(int a = 0; a < lattice.length; a++) {
				for(int b = 0; b < lattice[a].length; b++) {
					for(int c = 0; c < lattice[a][b].length; c++) {
						if(lattice[a][b][c] != null) {
							total2fold++;
							enerBefore = getS().calculateEnergyofFirstShell(lattice[a][b][c]);
							br1 = r.nextInt(4);
							while((br2=r.nextInt(4))!= br1) {;}
							cPos = lattice[a][b][c].center.getPosition();
							brPos = JVector.multiply(JVector.add(lattice[a][b][c].ligands[br1].getPosition(), lattice[a][b][c].ligands[br2].getPosition()), .5);
							translate = JVector.multiply(JVector.subtract(brPos, cPos), .5);
							lattice[a][b][c].translate(translate);
							enerAfter = getS().calculateEnergyofFirstShell(lattice[a][b][c]);
							if(enerAfter >= enerBefore) {
								e = -Math.abs((enerAfter-enerBefore))/ kb/T;
								e = Math.exp(e);
								if(e < r.nextDouble()) {
									lattice[a][b][c].translate(JVector.multiply(translate, -1));	
									twoFoldReset++;
								}
							}
						}
					}
				}
			}
		}
		System.out.println("Total number of 3-fold alignments: " + total3fold + "\t3-fold reset: " + threeFoldReset);
		System.out.println("Total number of 2-fold alignments: " + total2fold + "\t2-fold reset: " + twoFoldReset);
		return lattice;
	}
	public IdealTetrahedron[][][] monteCarloWalk(IdealTetrahedron[][][] lattice, JVector[] positions, int cycles, int idx) {
		//IdealTetrahedron[][][] lattice = makeTetragonalFCCLattice();
		getS().setLattice3(lattice);
		numMols = (int) (Math.pow(lattice.length/2, 3) * 4);
		getS().setLattice3(lattice);
		double enerBefore, enerAfter;
		int newChoice, oldChoice, twoOrThree;
		int total3fold = 0;
		int threeFoldReset = 0;
		int total2fold = 0;
		int twoFoldReset = 0;
		int distCutoff = 0;
		
		JVector move = new JVector();
		double e;
		int a, b, c;
		boolean newOrientation;
		
		//int walkTime = cycles / latticeCycles;
		for(int i = 0; i < cycles; i++) {
			/*if(cycles % walkTime == 0) {
				s.walk(walkLength, 1);
			}*/
			for(int j = 0; j < positions.length; j++) {
				a = (int) (Math.rint(positions[j].i));
				b = (int) (Math.rint(positions[j].j));
				c = (int) (Math.rint(positions[j].k));
				if(lattice[a][b][c] != null) {
					move.i = lattice[a][b][c].center.getPosition().i;
					move.j = lattice[a][b][c].center.getPosition().j;
					move.k = lattice[a][b][c].center.getPosition().k;
					newOrientation = checkDistances(lattice[a][b][c], distanceCutoff);
					total2fold++;
					if(newOrientation) {
						distCutoff++;
						enerBefore = getS().calculateEnergyofFirstShell(lattice[a][b][c]);
						oldChoice = lattice[a][b][c].orientation;
						twoOrThree = 0;
						switch(twoOrThree) {
						case 0:
							total2fold++;
							newChoice = r.nextInt(6);
							while(newChoice == oldChoice) { newChoice = r.nextInt(6); }
							lattice[a][b][c] = getOrientation(newChoice);
							lattice[a][b][c].translate(move);
							break;
						case 1:		
							total3fold++;
							if(idx == 1)
								alignSecondaryBromine(lattice, a, b, c, alignPrimaryBromine(lattice, a, b, c));
							else
								alignPrimaryBromine(lattice, a, b, c);
							break;
						}
						enerAfter = getS().calculateEnergyofFirstShell(lattice[a][b][c]);
						if(enerAfter >= enerBefore) {
							e = -Math.abs((enerAfter-enerBefore))/ kb/T;
							e = Math.exp(e);
							if(e < r.nextDouble()) {
								lattice[a][b][c] = getOrientation(oldChoice);
								lattice[a][b][c].translate(move);	
								if(twoOrThree == 1) { threeFoldReset++; }
								else { twoFoldReset++; }
									
							}
						}
					}
				}
			}
		}
		System.out.println("Total number of 2-fold tests: " + total2fold + "\tdistanceCutoff " + distCutoff + "\t2-fold reset: " + twoFoldReset);
		System.out.println("Total number of changed orientations: " + (distCutoff-twoFoldReset));
		return lattice;
	}
	public JVector getAlignedAxis(int orientation) {
		switch(orientation) {
		case 0:
		case 1:
			return JVector.x;

		case 2:
		case 3:
			return JVector.y;

		case 4:
		case 5:
			return JVector.z;
		}
		return null;
	}
	public JVector[] getOrthogonal110s(JVector twoFold) {
		Stack<Integer> s = new Stack<Integer>();
		int stackSize = 0;
		
		for(int i = 0; i < JVector.v110s.length; i++) {
			if(Math.rint(JVector.angle(twoFold, JVector.v110s[i])) == 90) { s.push(i); }
		}
		stackSize = s.size();
		JVector[] orthogonal110s = new JVector[stackSize];
		for(int i = 0; i < stackSize; i++) {
			orthogonal110s[i] = JVector.v110s[s.pop()];
		}
		return orthogonal110s;
	}
	public IdealTetrahedron[][][] monteCarloShift(IdealTetrahedron[][][] lattice, JVector[] positions, int cycles, int idx) {
		getS().setLattice3(lattice);
		numMols = (int) (Math.pow(lattice.length/2, 3) * 4);
		getS().setLattice3(lattice);
		double enerBefore, enerAfter;
		int total2fold = 0;
		int twoFoldReset = 0;
		
		JVector translate = new JVector();
		JVector cPos;
		JVector brPos;
		JVector twoFold;
		JVector[] orthogonal110s;
		JVector shift;
		int br1, br2;
		double e;
		int a, b, c;
		//int walkTime = cycles / latticeCycles;
		for(int i = 0; i < cycles; i++) {
			/*if(cycles % walkTime == 0) {
				s.walk(walkLength, 1);
			}*/
			for(int j = 0; j < positions.length; j++) {
				a = (int) (Math.rint(positions[j].i));
				b = (int) (Math.rint(positions[j].j));
				c = (int) (Math.rint(positions[j].k));
				if(lattice[a][b][c] != null) {
					total2fold++;
					enerBefore = getS().calculateEnergyofFirstShell(lattice[a][b][c]);
					twoFold = getAlignedAxis(lattice[a][b][c].orientation);
					if(twoFold == null) { 
						System.out.println("Two fold is null. Check your code.");
						continue; 
					}
					orthogonal110s = getOrthogonal110s(twoFold);
					shift = orthogonal110s[r.nextInt(orthogonal110s.length)];
					translate = JVector.multiply(shift, shift440);
					/*
					br1 = r.nextInt(4);
					while((br2=r.nextInt(4))!= br1) {;}
					cPos = lattice[a][b][c].center.position;
					brPos = JVector.multiply(JVector.add(lattice[a][b][c].ligands[br1].position, lattice[a][b][c].ligands[br2].position), .5);
					translate = JVector.multiply(JVector.subtract(brPos, cPos), .5);*/
					lattice[a][b][c].translate(translate);
					enerAfter = getS().calculateEnergyofFirstShell(lattice[a][b][c]);
					if(enerAfter >= enerBefore) {
						e = -Math.abs((enerAfter-enerBefore))/ kb/T;
						e = Math.exp(e);
						if(e < r.nextDouble()) {
							lattice[a][b][c].translate(JVector.multiply(translate, -1));	
							twoFoldReset++;
						}
					}
				}
			}
		}
		System.out.println("Total number of 2-fold shifts: " + total2fold + "\t2-fold reset: " + twoFoldReset);
		return lattice;
	}
	public void buildFCCLatticeWTetrahedron(int scattererType) {
		lattice = new JAtom[(int) (Math.pow(numUnits, 3) * 4*4)];
		int scattererIdx = 0;
		double cbrLen = 1.91/8.82;
		JAtom[] tetrahedron;
		JVector pos;
		JVector[] tetra = new JVector[4];
		tetra[0] = new JVector(cbrLen, cbrLen, cbrLen);
		tetra[1] = new JVector(-cbrLen, -cbrLen, cbrLen);
		tetra[2] = new JVector(-cbrLen, cbrLen, -cbrLen);
		tetra[3] = new JVector(cbrLen, -cbrLen, -cbrLen);
		for(int a = 0; a < numUnits*2; a++) {
			for(int b = 0; b < numUnits*2; b++) {
				for(int c = 0; c < numUnits*2; c++) {
					if((a + b + c) % 2 == 0) {
						pos = new JVector(a, b, c);
						pos = JVector.multiply(pos, .5);
						for(int j = 0; j < 4; j++) {
							lattice[scattererIdx] = new JAtom(scattererType, JVector.add(pos, tetra[j]));
							scattererIdx++;
						}
					}
				}
			}
		}
		System.out.println("num scatterers: " + scattererIdx);
	}

	public void makeAtomLatticeFromXYZ(File inpFile) {
		FileReader fr = null;

		Scanner s = null;
		
		int howManyLines = 0;
		
		try {
			fr = new FileReader(inpFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		double aDim = 1;
		s = new Scanner(fr);
		howManyLines = s.nextInt();

		s.nextLine();
		//s.nextLine();
		
		// define the lattice
		lattice = new JAtom[(howManyLines)];

		// make the lattice
		//s.nextLine();	// skip the CELL line
		s.nextLine();	// skip the SPGP line
		JAtom temp = null;
		for(int i = 0; i < lattice.length; i++)
		{
			try {
				temp = convertFromXYZLine(s.nextLine());
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			lattice[i] = (JAtom) temp.clone();
		}
	}
	public IdealTetrahedron[][][] dim1Todim3(IdealTetrahedron[] initLattice) {
		int i, j, k;
		int max=1;
		JVector pos;
		Scanner s = new Scanner(System.in);
		for(int a = 0; a < initLattice.length; a++) {
			pos = initLattice[a].getCenter().getPosition();
			pos = JVector.multiply(pos, 1./4.41);
			if(pos.i>max)
				max = (int)(Math.rint(pos.i));
			if(pos.j>max)
				max = (int)(Math.rint(pos.j));
			if(pos.k>max)
				max = (int)(Math.rint(pos.k));
		}
		max++;
		System.out.println("max: " + max);
		//s.next();
		IdealTetrahedron[][][] newLattice = new IdealTetrahedron[max][max][max];
		for(int a = 0; a < initLattice.length; a++) {
			pos = initLattice[a].getCenter().getPosition();
		//	System.out.print("oldPos: \t" + pos);
			pos = JVector.multiply(pos, 1./4.41);
			i = (int)Math.rint(pos.i);
			j = (int)Math.rint(pos.j);
			k = (int)Math.rint(pos.k);
			//System.out.println("\tnewPos:\t" + i + "," + j + "," + k);
			newLattice[i][j][k] = (IdealTetrahedron) initLattice[a].clone();
		}
		return newLattice;
	}
	public void readBoxes(File aFile, int whichBox) throws FileNotFoundException{
		FileInputStream fis = new FileInputStream(aFile);
		Scanner s = new Scanner(fis);
		int lines = s.nextInt();
		s.nextLine();
		int numBoxes = s.nextInt();
		s.nextLine();
		int tetraPerBox = s.nextInt();
		s.nextLine();

		IdealTetrahedron[][] boxes;
		
		String line = "";
		String[] splitLine;
		int idx = 0;
		JAtom[] atoms = new JAtom[5];
		boxes = new IdealTetrahedron[numBoxes][tetraPerBox];
		int tetraIdx = 0;
		int boxIdx = 0;
		int atomIdx = 0;
		int atomsPerTetrahedron = 5;
		while(s.hasNextLine()) {
			line = s.nextLine();
			splitLine = line.split("\t");
			
			if(splitLine.length == 3)
				// The line corresponds to the alignment vector
				continue;
			if(splitLine.length == 1)
				if(splitLine[0].compareTo("null") == 0)
					// The line is "null" because the algorithm did not find a tetrahedron at that location
					boxes[boxIdx][tetraIdx++] = null;
				else {
					// or the line is the new box index
					boxIdx = Integer.valueOf(splitLine[0]);
					tetraIdx = 0;
				}
			else if(splitLine.length == 4)
				// The line corresponds to an atom
				atoms[atomIdx++] = convertFromXYZLine(line);
			
			// check to see if I have 5 atoms yet
			if(atomIdx == atomsPerTetrahedron) {
				// make a new tetrahedron
				boxes[boxIdx][tetraIdx++] = new IdealTetrahedron(atoms);
				atomIdx = 0;
				atoms = null;
				atoms = new JAtom[5];
			}
				
		}
		System.out.println("Tetra per box: " + tetraPerBox);
		switch(whichBox) {
		case 1: boxes1 = new IdealTetrahedron[boxes.length][boxes[0].length];break;
		case 2: boxes2 = new IdealTetrahedron[boxes.length][boxes[0].length];break;
		case 3: boxes3 = new IdealTetrahedron[boxes.length][boxes[0].length];break;
		case 4: boxes4 = new IdealTetrahedron[boxes.length][boxes[0].length];break;
		case 5: boxes5 = new IdealTetrahedron[boxes.length][boxes[0].length];break;
		case 6: boxes6 = new IdealTetrahedron[boxes.length][boxes[0].length];break;
		case 7: boxes7 = new IdealTetrahedron[boxes.length][boxes[0].length];break;
		case FIRST_SHELLS: firstShells = new IdealTetrahedron[boxes.length][boxes[0].length];break;
		case SECOND_SHELLS: secondShells = new IdealTetrahedron[boxes.length][boxes[0].length];break;
		case STARS: stars = new IdealTetrahedron[boxes.length][boxes[0].length]; break;
		default: boxes2 = new IdealTetrahedron[boxes.length][boxes[0].length];break;
		}
		
		for(int i = 0; i < boxes.length; i++) {
			for(int j = 0; j < boxes[i].length; j++) {
				
				if(boxes[i][j] == null) 
					continue;
				
				switch(whichBox) {
				case 1: boxes1[i][j] = (IdealTetrahedron) boxes[i][j].clone(); break;
				case 2: boxes2[i][j] = (IdealTetrahedron) boxes[i][j].clone(); break;
				case 3: boxes3[i][j] = (IdealTetrahedron) boxes[i][j].clone(); break;
				case 4: boxes4[i][j] = (IdealTetrahedron) boxes[i][j].clone(); break;
				case 5: boxes5[i][j] = (IdealTetrahedron) boxes[i][j].clone(); break;
				case 6: boxes6[i][j] = (IdealTetrahedron) boxes[i][j].clone(); break;
				case 7: boxes7[i][j] = (IdealTetrahedron) boxes[i][j].clone(); break;
				case FIRST_SHELLS: firstShells[i][j]  = (IdealTetrahedron) boxes[i][j].clone(); break;
				case SECOND_SHELLS: secondShells[i][j]  = (IdealTetrahedron) boxes[i][j].clone(); break;
				case STARS: stars[i][j] = (IdealTetrahedron) boxes[i][j].clone(); break;
				default: boxes2[i][j] = (IdealTetrahedron) boxes[i][j].clone(); break;
				}
			}
		}
	}
	public IdealTetrahedron[] makeTetrahedronLatticeFromXYZ(File inpFile) {
		IdealTetrahedron[] lattice;
		FileReader fr = null;

		Scanner s = null;
		
		int howManyLines = 0;
		
		try {
			fr = new FileReader(inpFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		s = new Scanner(fr);
		howManyLines = s.nextInt();

		s.nextLine();
		//s.nextLine();
		double aDim = 1;	// read in the dimension of the box
		// define the lattice
		lattice = new IdealTetrahedron[(howManyLines-2)/5];

		// make the lattice
		//s.nextLine();	// skip the CELL line
		s.nextLine();	// skip the SPGP line
		boolean hasNextTetrahedron = true;
		for(int i = 0; i < lattice.length; i++)
		{
			JAtom C = null;
			hasNextTetrahedron = true;
			try {
				C = convertFromXYZLine(s.nextLine());
			} catch (Exception e1) {
				hasNextTetrahedron = false;
			}
			if(hasNextTetrahedron) {
				C.multiply(aDim);
				JAtom[] ligands = new JAtom[4];
				for(int j = 0; j < 4; j++)
				{
					try {
						ligands[j] = convertFromXYZLine(s.nextLine());
					} catch (Exception e) {
						e.printStackTrace();
					}
					ligands[j].multiply(aDim);
				}
				IdealTetrahedron toPlace = new IdealTetrahedron(C, ligands);
				
				lattice[i] = (IdealTetrahedron) toPlace.clone();
			}
		}
		
		return lattice;
	}
	
	public JAtom convertFromXYZLine(String xyzLine) {
		String[] temp = xyzLine.split("\t");
		
		if(temp[0].isEmpty())
			return null;

		int Z = Integer.valueOf(temp[0]);
		
		double x = Double.valueOf(temp[1]);
		
		double y = Double.valueOf(temp[2]);
		
		double z = Double.valueOf(temp[3]);

		return new JAtom(Z, new JVector(x, y, z));
	}
	public JAtom[] getLattice() { return lattice; }
	public void setRandom(Random r) { this.r = r; }

	public Simulate getS() {
		return s;
	}

	public void setS(Simulate s) {
		this.s = s;
	}
}
