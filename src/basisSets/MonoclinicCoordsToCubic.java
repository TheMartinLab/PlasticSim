package basisSets;

import jama.Matrix;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Random;
import java.util.Scanner;
import java.util.Stack;
import java.util.TreeMap;

import defaultPackage.DoubleLinkedListCBr4;
import defaultPackage.DoubleLinkedListVector;
import defaultPackage.IdealTetrahedron;
import defaultPackage.JAtom;
import defaultPackage.JVector;
import defaultPackage.Lattice;
import defaultPackage.LennardJonesLookup;
import defaultPackage.LennardJonesPotential;
import defaultPackage.TwoNodeCBr4;
import defaultPackage.TwoNodeVector;

public class MonoclinicCoordsToCubic {

	static double monoA = 21.43;
	static double monoB = 12.12;
	static double monoC = 21.02;
	protected int aUnits;
	protected double aSmall, aBig, aPlastic, aMono;
	protected JAtom[] carbons, bromines, box;
	protected IdealTetrahedron[] pseudoLattice, bigLattice, monoCoords, moreMonoCoords, cartesianCoords;
	protected IdealTetrahedron[][] pseudoLattices, miniLattices, fourCubes, fourBigCubes, lotsOLittleCubes, rotatedBigCubes, eightCubes,
		rotatedEightCubes, truncatedEightCubes;
	static double[][] m1 = makeTransformMatrix1();
	static double[][] m11 = makeTransformMatrix2();
	static JVector[] uniqueCenters_mono = getUniqueCenters_mono();
	static JVector[] uniqueCenters_cartesian = getUniqueCenters_cartesian();
	/**
	 * @return The eight unique positions in monoclinic coordinates
	 */
	public static JVector[] getUniqueCenters_mono() {
		JVector[] centers = new JVector[32];
		/*centers[0] = new JVector(.126, .316, .123);
		centers[1] = new JVector(.155, .719, .129);
		centers[2] = new JVector(.096, .032, .378);
		centers[3] = new JVector(.121, .559, .38 );
		centers[4] = new JVector(.155, .209, .629);
		centers[5] = new JVector(.126, .684, .623);
		centers[6] = new JVector(.121, .441, .88 );
		centers[7] = new JVector(.096, .968, .878);
		centers[8] = new JVector(.626, .816, .123);*/
		centers[0] = new JVector(0.345,0.709,0.871);
		centers[1] = new JVector(0.845,0.209,0.871);
		centers[2] = new JVector(0.655,0.709,0.629);
		centers[3] = new JVector(0.155,0.209,0.629);
		centers[4] = new JVector(0.655,0.291,0.129);
		centers[5] = new JVector(0.155,0.791,0.129);
		centers[6] = new JVector(0.345,0.291,0.371);
		centers[7] = new JVector(0.845,0.791,0.371);
		centers[8] = new JVector(0.626,0.184,0.623);
		centers[9] = new JVector(0.126,0.684,0.623);
		centers[10] = new JVector(0.374,0.184,0.877);
		centers[11] = new JVector(0.874,0.684,0.877);
		centers[12] = new JVector(0.374,0.816,0.377);
		centers[13] = new JVector(0.874,0.316,0.377);
		centers[14] = new JVector(0.626,0.816,0.123);
		centers[15] = new JVector(0.126,0.316,0.123);
		centers[16] = new JVector(0.379,0.941,0.62);
		centers[17] = new JVector(0.879,0.441,0.62);
		centers[18] = new JVector(0.621,0.941,0.88);
		centers[19] = new JVector(0.121,0.441,0.88);
		centers[20] = new JVector(0.621,0.059,0.38);
		centers[21] = new JVector(0.121,0.559,0.38);
		centers[22] = new JVector(0.379,0.059,0.12);
		centers[23] = new JVector(0.879,0.559,0.12);
		centers[24] = new JVector(0.596,0.468,0.878);
		centers[25] = new JVector(0.096,0.968,0.878);
		centers[26] = new JVector(0.404,0.468,0.622);
		centers[27] = new JVector(0.904,0.968,0.622);
		centers[28] = new JVector(0.404,0.532,0.122);
		centers[29] = new JVector(0.904,0.032,0.122);
		centers[30] = new JVector(0.596,0.532,0.378);
		centers[31] = new JVector(0.096,0.032,0.378);

		return centers;
	}
	public static JVector[] getUniqueCenters_cartesian() {
		JVector[] centers = new JVector[uniqueCenters_mono.length];
		for(int i = 0; i < uniqueCenters_mono.length; i++) {
			centers[i] = fractionalMonoToCartesianCubic(uniqueCenters_mono[i]);
			centers[i] = JVector.multiply(centers[i], 17.2*17.2);
		}
		return centers;
	}
	public void readFile2(File aFile) throws FileNotFoundException {
		FileInputStream fis = new FileInputStream(aFile);
		Scanner s = new Scanner(fis);
		
		int cIdx, brIdx, monoIdx;
		// figure out how many bromines and carbons I have

		carbons = new JAtom[32];
		bromines = new JAtom[32*4];
		monoCoords = new IdealTetrahedron[32];
		JAtom[] tempBr;
		cIdx = brIdx = monoIdx = 0;
		while(s.hasNextLine())
		{
			String[] next = s.nextLine().split("\t");
			int Z = Integer.valueOf(next[2]);
			
			double x = Double.valueOf(next[3]);
			double y = Double.valueOf(next[4]);
			double z = Double.valueOf(next[5]);
			
			if(Z == 6)
			{
				carbons[cIdx] = new JAtom(Z, new JVector(x, y, z));
				cIdx++;
			}
			
			if(Z == 35)
			{
				bromines[brIdx] = new JAtom(Z, new JVector(x, y, z));
				brIdx++;
			}
		}
		int newBromines = 0;
		for(int a = 0; a < bromines.length; a++) {
			for(int i = -1; i <= 1; i++) {
				for(int j = -1; j <= 1; j++) {
					for(int k = -1; k <= 1; k++) {
						newBromines++;
						
					}
				}
			}
		}
		tempBr = new JAtom[newBromines];
		newBromines = 0;
		for(int a = 0; a < bromines.length; a++) {
			for(int i = -1; i <= 1; i++) {
				for(int j = -1; j <= 1; j++) {
					for(int k = -1; k <= 1; k++) {
						tempBr[newBromines] = (JAtom) bromines[a].clone();
						tempBr[newBromines].setPosition(JVector.add(new JVector(i, j, k), tempBr[newBromines].getPosition()));
						newBromines++;
					}
				}
			}
		}
		JAtom[] ligands = new JAtom[4];
		double dist = 0;
		JAtom center;
		for(monoIdx = 0; monoIdx < carbons.length; monoIdx++) {
			DoubleLinkedListVector list = new DoubleLinkedListVector(5);
			center = (JAtom) carbons[monoIdx].clone();
			for(brIdx = 0; brIdx < tempBr.length; brIdx++) {
				dist = JVector.distance(center.getPosition(), tempBr[brIdx].getPosition());
				list.insert(new TwoNodeVector(dist, brIdx, null, null));
			}
			for(int i = 0; i < 4; i++) {
				ligands[i] = (JAtom) tempBr[list.removeHead().getKey()].clone();
			}
			monoCoords[monoIdx] = new IdealTetrahedron(center, (JAtom[]) ligands.clone());
			
		}
	}
	public void extendCrystal(int size) {
		int monoIdx = 0;
		for(int a = 0; a < monoCoords.length; a++) {
			for(int i = -size; i <= size; i++) {
				for(int j = -size; j <= size; j++) {
					for(int k = -size; k <= size; k++) {
						monoIdx++;
					}
				}
			}
		}
		
		moreMonoCoords = new IdealTetrahedron[monoIdx];
		monoIdx = 0;
		for(int a = 0; a < 32; a++) {
			for(int i = -size; i <= size; i++) {
				for(int j = -size; j <= size; j++) {
					for(int k = -size; k <= size; k++) {
						JVector pos = new JVector(i, j, k);
						moreMonoCoords[monoIdx] = (IdealTetrahedron) monoCoords[a].clone();
						moreMonoCoords[monoIdx].translate(pos);
						monoIdx++;
					}
				}
			}
		}
	}
	
	public void moreMonoCoordsToCartesianCubic() {
		cartesianCoords = new IdealTetrahedron[moreMonoCoords.length];
		for(int i = 0; i < moreMonoCoords.length; i++) {
			cartesianCoords[i] = (IdealTetrahedron) moreMonoCoords[i].clone();
			cartesianCoords[i].getCenter().setPosition(fractionalMonoToCartesianCubic(cartesianCoords[i].getCenter().getPosition()));
			cartesianCoords[i].getCenter().setPosition(JVector.multiply(cartesianCoords[i].getCenter().getPosition(), 17.2*17.2));
			for(int j = 0; j < 4; j++) {
				cartesianCoords[i].getLigands()[j].setPosition(fractionalMonoToCartesianCubic(cartesianCoords[i].getLigands()[j].getPosition()));
				cartesianCoords[i].getLigands()[j].setPosition(JVector.multiply(cartesianCoords[i].getLigands()[j].getPosition(), 17.2*17.2));
			}
			cartesianCoords[i].translate(JVector.multiply(cartesianCoords[i].getCenter().getPosition(), 18.5/17.4-1));
			                                    
		}
	}

	public static JVector fractionalMonoToCartesianCubic(JVector pos) {
		double xNew = JVector.dot(pos, new JVector(Math.cos(20.88*Math.PI/180), 0, 0)) ;
		double yNew = JVector.dot(pos, new JVector(0, 1, 0)) ;
		double zNew = JVector.dot(pos, new JVector(-Math.sin(20.88*Math.PI/180), 0, 1));
		// convert to cartesian coordinates
		xNew *= monoA;
		yNew *= monoB;
		zNew *= monoC;
		// convert to the new box coordinate system
		double[][] m2 = new double[4][1];
		
		m2[0][0] = xNew;
		m2[1][0] = yNew;
		m2[2][0] = zNew;
		m2[3][0] = 1;

		//TODO Check this new code. Might be error prone
		Matrix M1 = new Matrix(m1);
		Matrix M2 = new Matrix(m2);
		Matrix M11 = new Matrix(m11);
		
		Matrix MCoords = M1.times(M2);
		
		MCoords = M11.times(MCoords);
		
		double[][] newCoords = MCoords.getArray();
		
		return new JVector(newCoords[0][0], newCoords[1][0], newCoords[2][0]);
	}

	public static double[][] makeTransformMatrix1()
	{
		double[][] m1 = new double[4][4];
		m1[0][0] = 0.032964475;
		m1[0][1] = -0.032173565;
		m1[0][2] = 0.03338877;
		m1[0][3] = 0.311383956;
		m1[1][0] = -0.040505491;
		m1[1][1] = -8.12547E-05;
		m1[1][2] = 0.040664489;
		m1[1][3] = 0.39774325;
		m1[2][0] = 0.023218157;
		m1[2][1] = 0.047763466;
		m1[2][2] = 0.023517004;
		m1[2][3] = -0.625690357;
		m1[3][0] = 0;
		m1[3][1] = 0;
		m1[3][2] = 0;
		m1[3][3] = 1;
		return m1;
	}
	public static double[][] makeTransformMatrix2()
	{
		double[][] m1 = new double[4][4];
		m1[0][0] = 	0.056303224;
		m1[0][1] = -0.007155202;
		m1[0][2] = -0.010824335;
		m1[0][3] = 0;
		
		m1[1][0] = 0.009066961;
		m1[1][1] = 0.056528417;
		m1[1][2] = 0.006105919;
		m1[1][3] = 0;
		
		m1[2][0] = 0.010072661;
		m1[2][1] = -0.007247299;
		m1[2][2] = 0.055702074;
		m1[2][3] = 0;
		
		m1[3][0] = 0;
		m1[3][1] = 0;
		m1[3][2] = 0;
		m1[3][3] = 1;
		return m1;
	}
	public void idealizeTetrahedra(double newLen) {
		JAtom c;
		JAtom[] br;
		JVector oldPos, newPos;
		for(int i = 0; i < eightCubes.length; i++) {
			for(int j = 0; j < eightCubes[i].length; j++) {
				c = eightCubes[i][j].getCenter();
				br = eightCubes[i][j].getLigands();
				for(int a = 0; a < br.length; a++) {
					oldPos = JVector.subtract(c.getPosition(), br[a].getPosition());
					try {
						newPos = JVector.multiply(oldPos.unit(), newLen);
						br[a].setPosition(JVector.add(c.getPosition(), newPos));
					} catch (Exception e) {
					}
				}
			}
		}
	}
	public void makeEightCubes(int size) {
		int numPerCube = (int) (4*Math.pow(size, 3) + 6 * Math.pow(size, 2) + 3 * size + 1);
		System.out.println("numPerCube: " + numPerCube);
		eightCubes = new IdealTetrahedron[uniqueCenters_mono.length][63];
		for(int i = 0; i < eightCubes.length; i++) {
			eightCubes[i] = getCubeFromPosition(uniqueCenters_cartesian[i], size);
		}
	}
	public void rotateEightCubesx12() {
		rotatedEightCubes = new IdealTetrahedron[uniqueCenters_cartesian.length*12][eightCubes[0].length];
		JVector v111s = new JVector(1, 1, 1);
		double[] phi111 = {0, 120, 240};
		JVector[] v100s = {new JVector(0, 0, 0), new JVector(0, 0, 1), new JVector(0, 1, 0), new JVector(1, 0, 0)};
		double phi100 = 180;
		JVector origin = new JVector(aPlastic, aPlastic, aPlastic);
		int offset = 0;
		for(int a = 0; a < eightCubes.length; a++) {
			offset = 12*a;
			for(int b = 0; b < eightCubes[a].length; b++) {
				for(int i = 0; i < 3; i++) {
					rotatedEightCubes[offset+4*i][b] = (IdealTetrahedron) eightCubes[a][b].clone();
					rotatedEightCubes[offset+4*i][b].rotate(v111s, origin, phi111[i]);
					for(int j = 0; j < 4; j++) {
						rotatedEightCubes[offset+4*i+j][b] = (IdealTetrahedron) rotatedEightCubes[offset+4*i][b].clone();
						rotatedEightCubes[offset+4*i+j][b].rotate(v100s[j], origin, phi100);
					}
				}
			}
		}
	}
	public void zeroEightCubes() {
		DoubleLinkedListVector dllv = new DoubleLinkedListVector(rotatedEightCubes.length);
		JVector zero;
		double dist;
		for(int a = 0; a < rotatedEightCubes.length; a++) {
			// find the zero
			for(int b = 0; b < rotatedEightCubes[a].length; b++){ 
				dist = Math.abs(JVector.subtract(JVector.zero, rotatedEightCubes[a][b].getCenter().getPosition()).length());
				dllv.insert(new TwoNodeVector(dist, b, null, null));
			}
			// set the translation vector
			zero = JVector.multiply(rotatedEightCubes[a][dllv.removeHead().getKey()].getCenter().getPosition(), -1);
			// translate all molecules such that the zero is actually at zero
			for(int b = 0; b < rotatedEightCubes[a].length; b++){
				rotatedEightCubes[a][b].translate(zero);
			}
			dllv.clear();
		}
	}
	public void truncateEightx12Cubes(int size, int target) {
		// init the truncated eight cubes to be the same as the rotated eight cubes
		int endSize; 
		if(target%2 == 0) { endSize = (int) (4*Math.pow(target/2, 3) + 6 * Math.pow(target/2, 2) + 3*target/2 + 1); }
		else  { endSize = (int) (4*Math.pow((target+1)/2, 3)); }
		int numToRemove = 0;
		while(size > target) {
			size--;
			if(size % 2 == 0) {
				numToRemove += 2*Math.pow((size+1)/2+1, 2) + 2*((size+1)/2+1) + 1;
			}
			else {
				numToRemove += 2*Math.pow(size/2+1, 2) + 2*(size/2+1);
			}
		}
		
		IdealTetrahedron[][] temp = new IdealTetrahedron[rotatedEightCubes.length][rotatedEightCubes[0].length];
		truncatedEightCubes = new IdealTetrahedron[8*12][endSize];
		DoubleLinkedListVector dllv_Z = new DoubleLinkedListVector(rotatedEightCubes[0].length*2);
		DoubleLinkedListVector dllv_Y = new DoubleLinkedListVector(rotatedEightCubes[0].length*2);
		DoubleLinkedListVector dllv_X = new DoubleLinkedListVector(rotatedEightCubes[0].length*2);
		TwoNodeVector tnv;
		int idx;
		double val;
		for(int i = 0; i < rotatedEightCubes.length; i++) {
			idx = 0;
			for(int j = 0; j < rotatedEightCubes[i].length; j++) {
				temp[i][j] = (IdealTetrahedron) rotatedEightCubes[i][j].clone();
				val = JVector.dot(temp[i][j].getCenter().getPosition(), JVector.z);
				dllv_Z.insert(new TwoNodeVector(val, j, null, null));
				val = JVector.dot(temp[i][j].getCenter().getPosition(), JVector.y);
				dllv_Y.insert(new TwoNodeVector(val, j, null, null));
				val = JVector.dot(temp[i][j].getCenter().getPosition(), JVector.x);
				dllv_X.insert(new TwoNodeVector(val, j, null, null));
			}
			
			for(int j = 0; j < numToRemove; j++) {
				tnv = dllv_Z.removeLast();
				temp[i][tnv.getKey()] = null;
				tnv = dllv_Y.removeLast();
				temp[i][tnv.getKey()] = null;
				tnv = dllv_X.removeLast();
				temp[i][tnv.getKey()] = null;
			}
			for(int j = 0; j < temp[i].length; j++) {
				if(temp[i][j] == null) { continue; }
				truncatedEightCubes[i][idx] = (IdealTetrahedron) temp[i][j].clone();
				idx++;
			}
			dllv_Z.clear();
			dllv_Y.clear();
			dllv_X.clear();
		}
	}
	public IdealTetrahedron[][] truncateCubes(IdealTetrahedron[][] cubes, int size, int target) {
		// init the truncated eight cubes to be the same as the rotated eight cubes
		int endSize, numRemoved, startSize; 

		if(size%2 == 0) { startSize = (int) (4*Math.pow(size/2, 3) + 6 * Math.pow(size/2, 2) + 3*size/2 + 1); }
		else  { startSize = (int) (4*Math.pow((size+1)/2, 3)); }
		if(target%2 == 0) { endSize = (int) (4*Math.pow(target/2, 3) + 6 * Math.pow(target/2, 2) + 3*target/2 + 1); }
		else  { endSize = (int) (4*Math.pow((target+1)/2, 3)); }
		int numToRemove = 0;
		while(size > target) {
			
			if(size % 2 == 0) {
				numToRemove += Math.pow((size+1)/2+1, 2) + Math.pow((size+1)/2, 2);
			}
			else {
				numToRemove += 2*Math.pow((size+1)/2+1, 2);
			}
			size--;
		}
		numRemoved = 0;
		IdealTetrahedron[][] temp = new IdealTetrahedron[cubes.length][cubes[0].length];
		IdealTetrahedron[][] newCubes = new IdealTetrahedron[cubes.length][endSize];
		DoubleLinkedListVector dllv_Z = new DoubleLinkedListVector(cubes[0].length*2);
		DoubleLinkedListVector dllv_Y = new DoubleLinkedListVector(cubes[0].length*2);
		DoubleLinkedListVector dllv_X = new DoubleLinkedListVector(cubes[0].length*2);
		TwoNodeVector tnv;
		int idx;
		double val;
		for(int i = 0; i < cubes.length; i++) {
			idx = 0;
			for(int j = 0; j < cubes[i].length; j++) {
				temp[i][j] = (IdealTetrahedron) cubes[i][j].clone();
				val = JVector.dot(temp[i][j].getCenter().getPosition(), JVector.z);
				dllv_Z.insert(new TwoNodeVector(val, j, null, null));
				val = JVector.dot(temp[i][j].getCenter().getPosition(), JVector.y);
				dllv_Y.insert(new TwoNodeVector(val, j, null, null));
				val = JVector.dot(temp[i][j].getCenter().getPosition(), JVector.x);
				dllv_X.insert(new TwoNodeVector(val, j, null, null));
			}
			
			for(int j = 0; j < numToRemove; j++) {
				tnv = dllv_Z.removeLast();
				if(temp[i][tnv.getKey()] != null) {
					temp[i][tnv.getKey()] = null;
					numRemoved++;
					if(numRemoved == startSize-endSize) { break; }
				}
				
				tnv = dllv_Y.removeLast();
				if(temp[i][tnv.getKey()] != null) {
					temp[i][tnv.getKey()] = null;
					numRemoved++;
					if(numRemoved == startSize-endSize) { break; }
				}
				
				tnv = dllv_X.removeLast();
				if(temp[i][tnv.getKey()] != null) {
					temp[i][tnv.getKey()] = null;
					numRemoved++;
					if(numRemoved == startSize-endSize) { break; }
				}
			}
			for(int j = 0; j < temp[i].length; j++) {
				if(temp[i][j] == null) { continue; }
				newCubes[i][idx] = (IdealTetrahedron) temp[i][j].clone();
				idx++;
			}
			dllv_Z.clear();
			dllv_Y.clear();
			dllv_X.clear();
		}
		return newCubes;
	}
	public IdealTetrahedron[] getCubeFromPosition(JVector position, int size) {
		IdealTetrahedron[] cube = new IdealTetrahedron[ JVector.getFCCPositions(size).length];
		IdealTetrahedron zero;
		// find the zero point
		int numFinds = 0;
		double dist;
		JVector[] fccPos = JVector.getFCCPositions(size);
		DoubleLinkedListVector dllv = new DoubleLinkedListVector(10);
		for(int j = 0; j < fccPos.length; j++) {
			fccPos[j] = JVector.multiply(fccPos[j], aPlastic);
			for(int i = 0; i < cartesianCoords.length; i++) {
				dist = JVector.subtract(JVector.subtract(cartesianCoords[i].getCenter().getPosition(), fccPos[j]), position).length();
				dllv.insert(new TwoNodeVector(dist, i, null, null));
			}
			cube[j] = (IdealTetrahedron) cartesianCoords[dllv.removeHead().getKey()].clone();
			cube[j].translate(JVector.multiply(position, -1));
			dllv.clear();
		}
		return cube;
	}
	public JVector[] getFCCSites(int boxSize) {
		JVector[] sites = new JVector[(int) Math.pow(boxSize, 3)];
		int idx = 0;
		for(int i = 0; i < boxSize; i++) {
			for(int j = 0; j < boxSize; j++) {
				for(int k = 0; k < boxSize; k++) {
					if((i+j+k)%2 == 0) {
						sites[idx] = JVector.multiply(new JVector(i, j, k), aMono/4);
						idx++;
					}
				}
			}
		}
		return sites;
	}
	public void lookAtDisplacements(IdealTetrahedron[] box, int boxSize) {
		// get the six orientations
		IdealTetrahedron[] orientations = Lattice.getSix4barOrientations();
		DoubleLinkedListVector dllv = new DoubleLinkedListVector(6);
		DoubleLinkedListVector dllv2 = new DoubleLinkedListVector(6);
		JVector[] fcc = getFCCSites(boxSize);
		Double[] magnitudes = new Double[box.length];
		JVector[] displacements = new JVector[box.length];
		Double[] angles = new Double[box.length];
		JVector[] axes = JVector.v100s;
		JVector curPos;
		double angle = 0;
		TwoNodeVector curNode;
		double brDiff;
		double total = 0;
		JVector axis = null;
		Stack<Integer> s = new Stack<Integer>();
		System.out.println("shift\t\torientation\t\tangle\t\tbox[i].center.position\t\tcurPos");
		for(int i = 0; i < box.length; i++) {
			// figure out which fcc site the selected molecule is closest to
			curPos = JVector.multiply(JVector.multiply(box[i].getCenter().getPosition(), 4./aMono).roundInt(), aMono/4.);
			displacements[i] = JVector.subtract(box[i].getCenter().getPosition(), curPos);
			// figure out which 2-fold is aligned
			dllv.clear();
			for(int j = 0; j < orientations.length; j++) {
				orientations[j].moveTo(box[i].getCenter().getPosition());
				total = 0;
				for(int a = 0; a < 4; a++) {
					dllv2.clear();
					for(int b = 0; b < 4; b++) {
						brDiff = JVector.subtract(orientations[j].getLigands()[a].getPosition(), box[i].getLigands()[b].getPosition()).length();
						dllv2.insert(new TwoNodeVector(brDiff, b, null, null));
					}
					total += dllv2.removeHead().getValue();
				}
				dllv.insert(new TwoNodeVector(total, j, null, null));
			}
			int whichAxis = dllv.removeHead().getKey();
			switch(whichAxis)
			{
			case 0: axis = axes[0]; break;
			case 1: axis = axes[1]; break;
			case 2: axis = axes[2]; break;
			case 3: axis = axes[3]; break;
			case 4: axis = axes[4]; break;
			case 5: axis = axes[5]; break;
			}
			System.out.println(displacements[i]);
			// test which 110s are perpendicular to the chosen alignment vector
			for(int j = 0; j < JVector.v110s.length; j++) {
				brDiff = JVector.angle(axis, JVector.v110s[j]);
				if(Math.abs(brDiff-90) < 10) { s.push(j); }
			}
			
			System.out.println("\t\t100");
			for(int j = 0; j < JVector.v100s.length; j++) {
				if(j == whichAxis) { 
					System.out.print("*"); 
					System.out.println("\t\t" + JVector.v100s[j] + "\t" + JVector.angle(JVector.v100s[j], displacements[i]));
				}
			}
			System.out.println("\t\t110");
			for(int j = 0; j < JVector.v110s.length; j++) {
				if(s.contains(j)) { 
					System.out.print("*"); 
					System.out.println("\t\t" + JVector.v110s[j] + "\t" + JVector.angle(JVector.v110s[j], displacements[i]));
				}
			}
			System.out.println("\t" + box[i].getCenter().getPosition() + "\t" + curPos);
			s.clear();
		}
	}
	public void readFile(File aFile) throws FileNotFoundException
	{
		FileInputStream fis = new FileInputStream(aFile);
		Scanner s = new Scanner(fis);
		int num = 1331;

		int cIndex = 0;
		int brIndex = 0;
		int snIndex = 0;
		
		// figure out how many bromines and carbons I have
		while(s.hasNextLine())
		{
			String[] next = s.nextLine().split("\t");
			int Z = Integer.valueOf(next[2]);
			if(Z == 6)
				cIndex++;
			if(Z == 35)
				brIndex++;
			if(Z == 50)
				snIndex++;
		}
		carbons = new JAtom[num*32];
		bromines = new JAtom[num*4*128];
		box = new JAtom[snIndex];
		
		cIndex = 0;
		snIndex = 0;
		brIndex = 0;
		
		//double a = 21.43+ Math.sin(20.88*Math.PI/180)*21.02;
		//double b = 12.12;
		//double c = 21.02*(Math.cos(20.88*Math.PI/180));
		double a = 21.43;
		double b = 12.12;
		double c = 21.2;
		// reset the file input stream and the scanner
		fis = new FileInputStream(aFile);
		s = new Scanner(fis);
		//System.out.println("CELL\t" + a + "\t" + b + "\t" + c + "\t" + 90 + "\t" + 90 + "\t" + 90);
		//System.out.println("SPGP\tP1");
		
		double boxEdge = 17.4;
		int scalar = 1;
		while(s.hasNextLine())
		{
			if(snIndex == 10)
				System.out.print("");
			a = 21.43;
			b = 12.12;
			c = 21.02;
			
			String[] next = s.nextLine().split("\t");
			int Z = Integer.valueOf(next[2]);
			
			double x = Double.valueOf(next[3]);
			double y = Double.valueOf(next[4]);
			double z = Double.valueOf(next[5]);
			
			if(Z == 35) scalar = 4;
			else if(Z == 6) scalar = 1;
			JVector[] positions = new JVector[num*scalar];
			
			JVector pos = new JVector(x, y, z);

			int index = 0;

			int size = 5;
			for(int i = -size; i <= size; i++)
			{
				for(int j =-size; j <= size; j++)
				{
					for(int k = -size; k <= size; k++)
					{
						pos = JVector.add(new JVector(x, y, z), new JVector(i, j, k));
						// rotate the a axis
						pos = fractionalMonoToCartesianCubic(pos);
						
						if(Z == 50 && (i != 0 || j != 0 || k != 0))
							continue;
						
						JVector newPos = new JVector(pos.i, pos.j, pos.k);
						
						positions[index] = JVector.multiply(newPos, boxEdge*boxEdge);
						
						index++;
					}
				}
			}

			
			//JVector newPos2 = JVector.add(newPos, new JVector(a, 0, 0));
			//JVector newPos3 = JVector.add(newPos, new JVector(a, 0, c));
			//JVector newPos4 = JVector.add(newPos, new JVector(0, 0, c));
			
			// convert to new fractional coordinates
			
			a = 1;
			b = 1;
			c = 1;
			
			if(Z == 6)
			{
				for(int i = 0; i < positions.length; i++)
				{
					carbons[cIndex] = new JAtom(Z, positions[i]);
					cIndex++;
				}
			}
			
			if(Z == 35)
			{
				for(int i = 0; i < positions.length; i++)
				{			
					bromines[brIndex] = new JAtom(Z, positions[i]);
					brIndex++;
				}
			}
			if(Z == 50)
			{
				for(int i = 0; i < positions.length; i++)
				{
					if(positions[i]!= null && isInsideBox(positions[i], aBig))
					{
						box[snIndex] = new JAtom(Z, positions[i]);
						snIndex++;
					}
				}
			}
		}
		a = 21.43;
		b = 12.12;
		c = 21.02;
		
		makeTetrahedra();
	}

	public boolean isInsideBox(JVector position, double distance)
	{
		double aMin = -.027;
		double aMax = distance;
		double bMin = aMin;
		double bMax = aMax;
		double cMin = aMin;
		double cMax = aMax;
		double x = position.getI();
		double y = position.getJ();
		double z = position.getK();
		
		// test maximum positions
		if(x > aMax || y > bMax || z > cMax)
			return false;
		// test minimum positions
		if(x < aMin || y < bMin || z < cMin)
			return false;

		return true;
	}
	public void makeTetrahedra()
	{
 		int latticePos = 0;
		pseudoLattice = new IdealTetrahedron[8023*500];
		int ligandPos = 0;
		double distance = 0;
		// loop through the carbons
		for(int i = 0; i < carbons.length; i++)
		{
			
			JAtom[] ligands = new JAtom[4];
			ligandPos = 0;
			if(carbons[i] == null || carbons[i].getPosition() == null) { continue; }
			// loop through the bromines
			for(int j = 0; j < bromines.length; j++)
			{
				distance = 3;
				// calculate the distance from the bromine to the carbon
				if(bromines[j] == null || bromines[j].getPosition() == null) { continue; }
				
				distance = JVector.subtract(carbons[i].getPosition(), bromines[j].getPosition()).length();
				
				// add a bromine to the carbon tetrabromide tetrahedron
				if(distance < 2.5)
				{
					ligands[ligandPos] = (JAtom) bromines[j].clone();
					bromines[j] = null;
					ligandPos++;
					
					// if there are four ligands, move on to the next carbon center
					if(ligandPos == 4) { j = bromines.length; };
				}
			}
			if(ligandPos != 4) { continue; }
			// make a tetrahedron from the carbon and bromine positions
			IdealTetrahedron temp = new IdealTetrahedron(carbons[i], ligands);
			
			// calculate the new position for the molecule
			JVector newPos = JVector.multiply(temp.getCenter().getPosition(), aSmall/17.4);
			// subtract the old position from the new position
			newPos = JVector.subtract(newPos, temp.getCenter().getPosition());
			temp.translate(newPos);
			pseudoLattice[latticePos] = temp;
			latticePos++;
		}
	}
	public void rotatePseudoLattice()
	{
		pseudoLattices = new IdealTetrahedron[12][pseudoLattice.length];
		
		JVector[] translateBackIntoBox = new JVector[9];
		
		translateBackIntoBox[0] = new JVector(0, 1, 1);		// 110 - 100 - 180
		translateBackIntoBox[1] = new JVector(0, 0, 1);		// 110 - 100 - 270
		translateBackIntoBox[2] = new JVector(0, 0, 1);		// 110 - 010 - 90
		translateBackIntoBox[3] = new JVector(1, 0, 1);		// 110 - 010 - 180
		translateBackIntoBox[4] = new JVector(1, 1, 0);		// 110 - 001 - 180
		translateBackIntoBox[5] = new JVector(1, 0, 1);		// 101 - 001 - 180
		translateBackIntoBox[6] = new JVector(1, 0, 0);		// 101 - 001 - 270
		translateBackIntoBox[7] = new JVector(0, 1, 0);		// 101 - 001 - 270
		translateBackIntoBox[8] = new JVector(0, 1, 1);		// 011 - 100 - 180

		pseudoLattices = new IdealTetrahedron[12][pseudoLattice.length];
		// rotate around 111 by 0, 120 & 240 degrees
		for(int i = 0; i < pseudoLattice.length; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				IdealTetrahedron temp = (IdealTetrahedron) pseudoLattice[i].clone();
				temp.rotate(new JVector(1, 1, 1), new JVector(0, 0, 0), 120*j);
				pseudoLattices[j][i] = temp;
			}
		}
		// set the starting vector to be (110)
		int startingPoint = 0;
		JVector origin = new JVector(0, 0, 0);

		// rotate the 110 about the 100 by 180 & 270
		JVector axis = new JVector(1, 0, 0);
		for(int i = 0; i < pseudoLattice.length; i++)
		{
			for(int j = 2; j < 4; j++)
			{
				IdealTetrahedron temp = (IdealTetrahedron) pseudoLattices[startingPoint][i].clone();
				temp.rotate(axis, origin, 90*j);
				temp.translate(JVector.multiply(translateBackIntoBox[j-2], aBig));
				pseudoLattices[j+1][i] = temp;
			}
		}
		// rotate the 110 about the 010 by 90 & 180
		axis = new JVector(0, 1, 0);
		for(int i = 0; i < pseudoLattice.length; i++)
		{
			for(int j = 1; j < 3; j++)
			{
				IdealTetrahedron temp = (IdealTetrahedron) pseudoLattices[startingPoint][i].clone();
				temp.rotate(axis, origin, 90*j);
				temp.translate(JVector.multiply(translateBackIntoBox[j+1], aBig));
				pseudoLattices[j+4][i] = temp;

			}
		}
		// rotate the 110 about the 001 by 180
		axis = new JVector(0, 0, 1);
		for(int i = 0; i < pseudoLattice.length; i++)
		{
			IdealTetrahedron temp = (IdealTetrahedron) pseudoLattices[startingPoint][i].clone();
			temp.rotate(axis, origin, 180);
			temp.translate(JVector.multiply(translateBackIntoBox[4], aBig));
			pseudoLattices[7][i] = temp;

		}
		// set the starting vector to be (101)
		startingPoint = 1;
		// rotate the 101 about the 010 by 180, 270
		axis = new JVector(0, 1, 0);
		for(int i = 0; i < pseudoLattice.length; i++)
		{
			for(int j = 2; j < 4; j++)
			{
				IdealTetrahedron temp = (IdealTetrahedron) pseudoLattices[startingPoint][i].clone();
				temp.rotate(axis, origin, 90*j);
				temp.translate(JVector.multiply(translateBackIntoBox[j+3], aBig));
				pseudoLattices[j+6][i] = temp;
			}
		}
		// rotate the 101 about the 001 by 270
		axis = new JVector(0, 0, 1);
		for(int i = 0; i < pseudoLattice.length; i++)
		{
			IdealTetrahedron temp = (IdealTetrahedron) pseudoLattices[startingPoint][i].clone();
			temp.rotate(axis, origin, 270);
			temp.translate(JVector.multiply(translateBackIntoBox[7], aBig));
			pseudoLattices[10][i] = temp;

		}
		// set the starting vector to be (011)
		startingPoint = 2;
		// rotate the 011 about the 100 by 180
		axis = new JVector(1, 0, 0);
		for(int i = 0; i < pseudoLattice.length; i++)
		{
			IdealTetrahedron temp = (IdealTetrahedron) pseudoLattices[startingPoint][i].clone();
			temp.rotate(axis, origin, 180);
			temp.translate(JVector.multiply(translateBackIntoBox[8], aBig));
			pseudoLattices[11][i] = temp;
		}
	}
	
	public void truncatePseudoLattices()
	{
		IdealTetrahedron[][] newPseudoLattices = new IdealTetrahedron[12][32];
		
		for(int i = 0; i < pseudoLattices.length; i++)
		{
			int index = 0;
			for(int j = 0; j < pseudoLattices[i].length; j++)
			{
				if(pseudoLattices[i][j] != null && isInsideBox(pseudoLattices[i][j].getCenter().getPosition(), aSmall))
				{
					newPseudoLattices[i][index] = pseudoLattices[i][j];
					index++;
				}
			}
		}
		pseudoLattices = newPseudoLattices;
	}
	
	public void cutUpPseudoLattices()
	{
		JVector[] positions = new JVector[8];
		positions[0] = new JVector(0, 0, 0);
		positions[1] = new JVector(1, 0, 0);
		positions[2] = new JVector(0, 1, 0);
		positions[3] = new JVector(0, 0, 1);
		positions[4] = new JVector(1, 1, 0);
		positions[5] = new JVector(1, 0, 1);
		positions[6] = new JVector(0, 1, 1);
		positions[7] = new JVector(1, 1, 1);
		
		JVector offset = new JVector(2, 2, 2);
		for(int i = 0; i < positions.length; i++)
		{
			positions[i] = JVector.multiply(positions[i], aSmall/2);
			positions[i] = JVector.add(positions[i], offset);
		}
		
		//divide the 12 2x2x2 pseudo lattice cubes into 96 1x1x1 cubes
		
		// find the "origins" of the 8 1x1x1 cubes
		
		miniLattices = new IdealTetrahedron[96][4];

		int miniLatticeIndex = 0;
		
		// loop through the pseudoLattices
		for(int latticeIndex = 0; latticeIndex < pseudoLattices.length; latticeIndex++)
		{
			// get the eight centers
			for(int i = 0; i < positions.length; i++)
			{
				DoubleLinkedListCBr4 findCenter = new DoubleLinkedListCBr4(3);
				// find the center
				for(int j = 0; j < pseudoLattices[latticeIndex].length; j++)
				{
					if(pseudoLattices[latticeIndex][j] != null)
					{
						double howClose = Math.abs(JVector.subtract(pseudoLattices[latticeIndex][j].getCenter().getPosition(), positions[i]).length());
						findCenter.insert(new TwoNodeCBr4(pseudoLattices[latticeIndex][j], howClose, null, null));
					}

				}
				// get the center
				IdealTetrahedron center = findCenter.removeHead().getValue();
				
				// set the center value to null
				for(int j = 0; j < pseudoLattices[latticeIndex].length; j++)
				{
					if(pseudoLattices[latticeIndex][j] != null && pseudoLattices[latticeIndex][j].getID() == center.getID())
					{
						pseudoLattices[latticeIndex][j] = null;
						j = pseudoLattices[latticeIndex].length;
					}
				}
				
				// find the three closest molecules
				DoubleLinkedListCBr4 findClosest = new DoubleLinkedListCBr4(3);
				
				for(int j = 0; j < pseudoLattices[latticeIndex].length; j++)
				{
					if(pseudoLattices[latticeIndex][j] != null)
					{
						double howClose = Math.abs(JVector.subtract(pseudoLattices[latticeIndex][j].getCenter().getPosition(), center.getCenter().getPosition()).length());
						findClosest.insert(new TwoNodeCBr4(pseudoLattices[latticeIndex][j], howClose, null, null));
					}
				}
				miniLattices[miniLatticeIndex][0] = center;
				miniLattices[miniLatticeIndex][1] = findClosest.removeHead().getValue();
				miniLattices[miniLatticeIndex][2] = findClosest.removeHead().getValue();
				miniLattices[miniLatticeIndex][3] = findClosest.removeHead().getValue();
				
				// set the three closest values to null
				for(int j = 0; j < pseudoLattices[latticeIndex].length; j++)
				{
					if(pseudoLattices[latticeIndex][j] != null &&
							(miniLattices[miniLatticeIndex][1].getID() == pseudoLattices[latticeIndex][j].getID() ||
									miniLattices[miniLatticeIndex][2].getID() == pseudoLattices[latticeIndex][j].getID() ||
									miniLattices[miniLatticeIndex][3].getID() == pseudoLattices[latticeIndex][j].getID()))
					{
						pseudoLattices[latticeIndex][j] = null;
					}
				}
				
				// shift the mini lattice back into the correct size box
				JVector shift = JVector.subtract(positions[i], offset);
				for(int j = 0; j < miniLattices[miniLatticeIndex].length; j++)
					miniLattices[miniLatticeIndex][j].translate(JVector.multiply(shift, -1));
				miniLatticeIndex++;
			}
		}
		
	}
	public void getFourSmallCubes()
	{
		JVector[] positions = new JVector[8];
		positions[0] = new JVector(0, 0, 0);
		positions[1] = new JVector(1, 0, 0);
		positions[2] = new JVector(0, 1, 0);
		positions[3] = new JVector(0, 0, 1);
		positions[4] = new JVector(1, 1, 0);
		positions[5] = new JVector(1, 0, 1);
		positions[6] = new JVector(0, 1, 1);
		positions[7] = new JVector(1, 1, 1);
		
		JVector offset = new JVector(2, 2, 2);
		for(int i = 0; i < positions.length; i++)
		{
			positions[i] = JVector.multiply(positions[i], aSmall/2);
			positions[i] = JVector.add(positions[i], offset);
		}
		
		//divide the 12 2x2x2 pseudo lattice cubes into 96 1x1x1 cubes
		
		// find the "origins" of the 8 1x1x1 cubes
		
		miniLattices = new IdealTetrahedron[96][4];
	
		int miniLatticeIndex = 0;
		
		// get the eight centers
		for(int i = 0; i < positions.length; i++)
		{
			DoubleLinkedListCBr4 findCenter = new DoubleLinkedListCBr4(3);
			// find the center
			for(int j = 0; j < pseudoLattice.length; j++)
			{
				if(pseudoLattice[j] != null)
				{
					double howClose = Math.abs(JVector.subtract(pseudoLattice[j].getCenter().getPosition(), positions[i]).length());
					findCenter.insert(new TwoNodeCBr4(pseudoLattice[j], howClose, null, null));
				}

			}
			// get the center
			IdealTetrahedron center = findCenter.removeHead().getValue();
			
			// set the center value to null
			for(int j = 0; j < pseudoLattice.length; j++)
			{
				if(pseudoLattice[j] != null && pseudoLattice[j].getID() == center.getID())
				{
					pseudoLattice[j] = null;
					j = pseudoLattice.length;
				}
			}
			
			// find the three closest molecules
			DoubleLinkedListCBr4 findClosest = new DoubleLinkedListCBr4(3);
			
			for(int j = 0; j < pseudoLattice.length; j++)
			{
				if(pseudoLattice[j] != null)
				{
					double howClose = Math.abs(JVector.subtract(pseudoLattice[j].getCenter().getPosition(), center.getCenter().getPosition()).length());
					findClosest.insert(new TwoNodeCBr4(pseudoLattice[j], howClose, null, null));
				}
			}
			miniLattices[miniLatticeIndex][0] = center;
			miniLattices[miniLatticeIndex][1] = findClosest.removeHead().getValue();
			miniLattices[miniLatticeIndex][2] = findClosest.removeHead().getValue();
			miniLattices[miniLatticeIndex][3] = findClosest.removeHead().getValue();
			
			// set the three closest values to null
			for(int j = 0; j < pseudoLattice.length; j++)
			{
				if(pseudoLattice[j] != null &&
						(miniLattices[miniLatticeIndex][1].getID() == pseudoLattice[j].getID() ||
								miniLattices[miniLatticeIndex][2].getID() == pseudoLattice[j].getID() ||
								miniLattices[miniLatticeIndex][3].getID() == pseudoLattice[j].getID()))
				{
					pseudoLattice[j] = null;
				}
			}
			
			// shift the mini lattice back into the correct size box
			JVector shift = JVector.subtract(positions[i], offset);
			for(int j = 0; j < miniLattices[miniLatticeIndex].length; j++)
				miniLattices[miniLatticeIndex][j].translate(JVector.multiply(shift, -1));
			miniLatticeIndex++;
		}
		
		fourCubes = new IdealTetrahedron[4][4];
		
		fourCubes[0] = miniLattices[0];
		fourCubes[1] = miniLattices[4];
		fourCubes[2] = miniLattices[2];
		fourCubes[3] = miniLattices[7];
		
	}

	public void makeFourBigCubes()
	{
		IdealTetrahedron aa = fourCubes[0][0];
		IdealTetrahedron ab = fourCubes[0][2];
		IdealTetrahedron ac = fourCubes[0][1];
		IdealTetrahedron ad = fourCubes[0][3];
		
		IdealTetrahedron ba = fourCubes[1][0];
		IdealTetrahedron bb = fourCubes[1][3];
		IdealTetrahedron bc = fourCubes[1][1];
		IdealTetrahedron bd = fourCubes[1][2];
		
		IdealTetrahedron ca = fourCubes[2][0];
		IdealTetrahedron cb = fourCubes[2][3];
		IdealTetrahedron cc = fourCubes[2][2];
		IdealTetrahedron cd = fourCubes[2][1];
		
		IdealTetrahedron da = fourCubes[3][0];
		IdealTetrahedron db = fourCubes[3][1];
		IdealTetrahedron dc = fourCubes[3][2];
		IdealTetrahedron dd = fourCubes[3][3];
		
		fourBigCubes = new IdealTetrahedron[4][14];
		
		// box1
		// 0
		fourBigCubes[0][0] = (IdealTetrahedron) aa.clone();
		fourBigCubes[0][1] = (IdealTetrahedron) ca.clone();
		fourBigCubes[0][2] = (IdealTetrahedron) ab.clone();
		fourBigCubes[0][3] = (IdealTetrahedron) ca.clone();
		fourBigCubes[0][4] = (IdealTetrahedron) ba.clone();
		// 1/2
		fourBigCubes[0][5] = (IdealTetrahedron) ac.clone();
		fourBigCubes[0][6] = (IdealTetrahedron) ad.clone();
		fourBigCubes[0][7] = (IdealTetrahedron) cd.clone();
		fourBigCubes[0][8] = (IdealTetrahedron) cc.clone();
		// 1
		fourBigCubes[0][9] = (IdealTetrahedron) ca.clone();
		fourBigCubes[0][10] = (IdealTetrahedron) ba.clone();
		fourBigCubes[0][11] = (IdealTetrahedron) cb.clone();
		fourBigCubes[0][12] = (IdealTetrahedron) ba.clone();
		fourBigCubes[0][13] = (IdealTetrahedron) da.clone();
		
		// box2
		// 0
		fourBigCubes[1][0] = (IdealTetrahedron) ba.clone();
		fourBigCubes[1][1] = (IdealTetrahedron) da.clone();
		fourBigCubes[1][2] = (IdealTetrahedron) bb.clone();
		fourBigCubes[1][3] = (IdealTetrahedron) da.clone();
		fourBigCubes[1][4] = (IdealTetrahedron) aa.clone();
		// 1/2
		fourBigCubes[1][5] = (IdealTetrahedron) bc.clone();
		fourBigCubes[1][6] = (IdealTetrahedron) bd.clone();
		fourBigCubes[1][7] = (IdealTetrahedron) dd.clone();
		fourBigCubes[1][8] = (IdealTetrahedron) dc.clone();
		// 1
		fourBigCubes[1][9] = (IdealTetrahedron) da.clone();
		fourBigCubes[1][10] = (IdealTetrahedron) aa.clone();
		fourBigCubes[1][11] = (IdealTetrahedron) db.clone();
		fourBigCubes[1][12] = (IdealTetrahedron) aa.clone();
		fourBigCubes[1][13] = (IdealTetrahedron) ca.clone();
		
		// box3
		// 0
		fourBigCubes[2][0] = (IdealTetrahedron) ca.clone();
		fourBigCubes[2][1] = (IdealTetrahedron) ba.clone();
		fourBigCubes[2][2] = (IdealTetrahedron) cb.clone();
		fourBigCubes[2][3] = (IdealTetrahedron) ba.clone();
		fourBigCubes[2][4] = (IdealTetrahedron) da.clone();
		// 1/2
		fourBigCubes[2][5] = (IdealTetrahedron) cc.clone();
		fourBigCubes[2][6] = (IdealTetrahedron) cd.clone();
		fourBigCubes[2][7] = (IdealTetrahedron) bd.clone();
		fourBigCubes[2][8] = (IdealTetrahedron) bc.clone();
		// 1
		fourBigCubes[2][9] = (IdealTetrahedron) ba.clone();
		fourBigCubes[2][10] = (IdealTetrahedron) da.clone();
		fourBigCubes[2][11] = (IdealTetrahedron) bb.clone();
		fourBigCubes[2][12] = (IdealTetrahedron) da.clone();
		fourBigCubes[2][13] = (IdealTetrahedron) aa.clone();
		
		// box4
		// 0
		fourBigCubes[3][0] = (IdealTetrahedron) da.clone();
		fourBigCubes[3][1] = (IdealTetrahedron) aa.clone();
		fourBigCubes[3][2] = (IdealTetrahedron) db.clone();
		fourBigCubes[3][3] = (IdealTetrahedron) aa.clone();
		fourBigCubes[3][4] = (IdealTetrahedron) ca.clone();
		// 1/2
		fourBigCubes[3][5] = (IdealTetrahedron) dc.clone();
		fourBigCubes[3][6] = (IdealTetrahedron) dd.clone();
		fourBigCubes[3][7] = (IdealTetrahedron) ad.clone();
		fourBigCubes[3][8] = (IdealTetrahedron) ac.clone();
		// 1
		fourBigCubes[3][9] = (IdealTetrahedron) aa.clone();
		fourBigCubes[3][10] = (IdealTetrahedron) ca.clone();
		fourBigCubes[3][11] = (IdealTetrahedron) ab.clone();
		fourBigCubes[3][12] = (IdealTetrahedron) ca.clone();
		fourBigCubes[3][13] = (IdealTetrahedron) ba.clone();
		
		for(int i = 0; i < fourBigCubes.length; i++)
		{
			fourBigCubes[i][0].translate(new JVector(0, 0, 0));	//000
			fourBigCubes[i][1].translate(new JVector(aSmall/2, 0, 0));	// 100
			fourBigCubes[i][2].translate(new JVector(0, 0, 0));	// .5 .5 0
			fourBigCubes[i][3].translate(new JVector(0, aSmall/2, 0));	// 010
			fourBigCubes[i][4].translate(new JVector(aSmall/2, aSmall/2, 0));	//110
			
			fourBigCubes[i][5].translate(new JVector(0, 0, 0));	//.5 0 .5
			fourBigCubes[i][6].translate(new JVector(0, 0, 0));	// 0 .5 .5
			fourBigCubes[i][7].translate(new JVector(aSmall/2, 0, 0));	//1 .5 .5 
			fourBigCubes[i][8].translate(new JVector(0, aSmall/2, 0));	// .5 1 .5
			
			fourBigCubes[i][9].translate(new JVector(0, 0, aSmall/2));	// 001
			fourBigCubes[i][10].translate(new JVector(aSmall/2, 0, aSmall/2));	//101
			fourBigCubes[i][11].translate(new JVector(0, 0, aSmall/2));	//.5 .5 1
			fourBigCubes[i][12].translate(new JVector(0, aSmall/2, aSmall/2));	//011
			fourBigCubes[i][13].translate(new JVector(aSmall/2, aSmall/2, aSmall/2));	//111
		}
	}
	
	public void rotateBigCubes()
	{
		double a = aSmall*3/8;
		JVector origin = new JVector(a, a, a);
		
		int cubeIndex = 0;
		rotatedBigCubes = new IdealTetrahedron[48][fourBigCubes[0].length];
		// loop through all four cubes
		for(int i = 0; i < 4; i++)
		{
			// rotate the 100 to the 010 and 001 by rotating 120 and 240 degrees around the 111
			for(int j = 0; j < 3; j++)
			{
				// make the basis sets
				for(int k = 0; k < fourBigCubes[i].length; k++)
				{
					rotatedBigCubes[cubeIndex][k] = (IdealTetrahedron) fourBigCubes[i][k].clone();
					rotatedBigCubes[cubeIndex][k].rotate(new JVector(1, 1, 1), new JVector(0, 0, 0), 120*j);
				}
				cubeIndex++;
			}
			// rotate the 110 about the 100 by 180
			for(int k = 0; k < rotatedBigCubes[cubeIndex].length; k++)
			{
				rotatedBigCubes[cubeIndex][k] = (IdealTetrahedron) rotatedBigCubes[0][k].clone();
				rotatedBigCubes[cubeIndex][k].rotate(new JVector(1, 0, 0), origin, 180);
			}
			cubeIndex++;
			
			// rotate the 110 about the 100 by 270
			for(int k = 0; k < rotatedBigCubes[cubeIndex].length; k++)
			{
				rotatedBigCubes[cubeIndex][k] = (IdealTetrahedron) rotatedBigCubes[0][k].clone();
				rotatedBigCubes[cubeIndex][k].rotate(new JVector(1, 0, 0), origin, 270);
			}
			cubeIndex++;
			
			// rotate the 110 about the 010 by 90
			for(int k = 0; k < rotatedBigCubes[0].length; k++)
			{
				rotatedBigCubes[cubeIndex][k] = (IdealTetrahedron) rotatedBigCubes[0][k].clone();
				rotatedBigCubes[cubeIndex][k].rotate(new JVector(0, 1, 0), origin, 90);
			}
			cubeIndex++;
			
			// rotate the 110 about the 010 by 180
			for(int k = 0; k < rotatedBigCubes[0].length; k++)
			{
				rotatedBigCubes[cubeIndex][k] = (IdealTetrahedron) rotatedBigCubes[0][k].clone();
				rotatedBigCubes[cubeIndex][k].rotate(new JVector(0, 1, 0), origin, 180);
			}
			cubeIndex++;
			
			// rotate the 110 about the 001 by 180
			for(int k = 0; k < rotatedBigCubes[0].length; k++)
			{
				rotatedBigCubes[cubeIndex][k] = (IdealTetrahedron) rotatedBigCubes[0][k].clone();
				rotatedBigCubes[cubeIndex][k].rotate(new JVector(0, 0, 1), origin, 180);
			}
			cubeIndex++;
			
			// rotate the 101 about the 010 by 180, 270
			for(int k = 0; k < rotatedBigCubes[0].length; k++)
			{
				rotatedBigCubes[cubeIndex][k] = (IdealTetrahedron) rotatedBigCubes[1][k].clone();
				rotatedBigCubes[cubeIndex][k].rotate(new JVector(0, 1, 0), origin, 180);
			}
			cubeIndex++;
			
			// rotate the 101 about the 010 by 180, 270
			for(int k = 0; k < rotatedBigCubes[0].length; k++)
			{
				rotatedBigCubes[cubeIndex][k] = (IdealTetrahedron) rotatedBigCubes[1][k].clone();
				rotatedBigCubes[cubeIndex][k].rotate(new JVector(0, 1, 0), origin, 270);
			}
			cubeIndex++;
			
			// rotate the 101 about the 001 by 270
			for(int k = 0; k < rotatedBigCubes[0].length; k++)
			{
				rotatedBigCubes[cubeIndex][k] = (IdealTetrahedron) rotatedBigCubes[1][k].clone();
				rotatedBigCubes[cubeIndex][k].rotate(new JVector(0, 0, 1), origin, 270);
			}
			cubeIndex++;
			// rotate the 011 about the 100 by 180
			for(int k = 0; k < rotatedBigCubes[0].length; k++)
			{
				rotatedBigCubes[cubeIndex][k] = (IdealTetrahedron) rotatedBigCubes[2][k].clone();
				rotatedBigCubes[cubeIndex][k].rotate(new JVector(1, 0, 0), origin, 180);
			}
			cubeIndex++;
		}
	}
	
	public void truncateBigCubes()
	{
		lotsOLittleCubes = new IdealTetrahedron[48][4];
		JVector origin = new JVector(2, 2, 2);
		for(int i = 0; i < rotatedBigCubes.length; i++)
		{
			DoubleLinkedListCBr4 four = new DoubleLinkedListCBr4(6);
			for(int j = 0; j < rotatedBigCubes[i].length; j++)
			{
				double distance = JVector.subtract(rotatedBigCubes[i][j].getCenter().getPosition(), origin).length();
				four.insert(new TwoNodeCBr4(rotatedBigCubes[i][j], distance, null, null));
			}
			
			lotsOLittleCubes[i][0] = four.removeHead().getValue();
			lotsOLittleCubes[i][1] = four.removeHead().getValue();
			lotsOLittleCubes[i][2] = four.removeHead().getValue();
			lotsOLittleCubes[i][3] = four.removeHead().getValue();
		}
	}
	public void makeBigBox(int a, int b, int c)
	{
		int length = lotsOLittleCubes[0].length ;
		
		bigLattice = new IdealTetrahedron[length * a * b * c];
		
		//int whichRotation = 0;
		int index = 0;
		
		for(int i = 0; i < a; i++)
		{
			for(int j = 0; j < b; j++)
			{
				for(int k = 0; k < c; k++)
				{
					
					JVector translate = JVector.multiply(new JVector(i, j, k), aSmall/2);
					Random r = new Random();
					int whichLattice = r.nextInt(48);
					for(int m = 0; m < lotsOLittleCubes[0].length; m++)
					{
						IdealTetrahedron temp = (IdealTetrahedron) lotsOLittleCubes[whichLattice][m].clone();
						temp.translate(translate);
						bigLattice[index] = temp;
						index++;
					}
				}
			}
		}
	}
	
	public void makeRandomBigBox(int a, int b, int c)
	{
		int length = miniLattices[0].length;
		bigLattice = new IdealTetrahedron[length * a * b * c];

		int index = 0;
		
		for(int i = 0; i < a; i++)
		{
			for(int j = 0; j < b; j++)
			{
				for(int k = 0; k < c; k++)
				{
					
					JVector translate = JVector.multiply(new JVector(i, j, k), aSmall/2);
					Random r = new Random();
					int whichLattice = r.nextInt(12*8);
					for(int m = 0; m < miniLattices[0].length; m++)
					{
						IdealTetrahedron temp = (IdealTetrahedron) miniLattices[whichLattice][m].clone();
						temp.translate(translate);
						bigLattice[index] = temp;
						index++;
					}
				}
			}
		}
		
	}
	
	public void printPseudoLatticeToConsole()
	{
		/*
		double a = 21.43;
		double b = 12.12;
		double c = 21.02;
		*/
		double a = aBig/2;
		double b = a;
		double c = a;
		
		int aUnits = 1;
		int bUnits = 1;
		int cUnits = 1;
		System.out.println("CELL\t" + a + "\t" + b + "\t" + c + "\t" + 90 + "\t" + 90 + "\t" + 90);
		System.out.println("SPGP\tP1");				

		for(int i = 0; i < pseudoLattice.length; i++)
		{
			IdealTetrahedron temp = pseudoLattice[i];
			if(temp == null)
				continue;
			try
			{
				System.out.print(temp.toStringForAtoms(a, b, c, i, aUnits, bUnits, cUnits));
			}
			catch(NullPointerException npe) {}
			
		}
	}
	
	public void printBigBoxToConsole()
	{
		/*
		double a = 21.43;
		double b = 12.12;
		double c = 21.02;
		*/
		double a = 1;
		double b = a;
		double c = a;
		
		int aUnits = 1;
		int bUnits = 1;
		int cUnits = 1;
		System.out.println("CELL\t" + a + "\t" + b + "\t" + c + "\t" + 90 + "\t" + 90 + "\t" + 90);
		System.out.println("SPGP\tP1");				

		for(int i = 0; i < bigLattice.length; i++)
		{
			IdealTetrahedron temp = bigLattice[i];
			if(temp == null)
				continue;
			try
			{
				System.out.print(temp.toStringForAtoms(a, b, c, i, aUnits, bUnits, cUnits));
			}
			catch(NullPointerException npe) {}
			
		}
	}
	
	public String makeAtomsHeader(double a)
	{
		String header = "";
		
		//double a = aUnits * aSmall/2;
		
		header += "CELL " + a + " " + a + " " + a + " " + 90 + " " + 90 + " " + 90 + "\n";
		header += "SPGP P1\n";
		
		return header;
	}
	
	public void printToFile(File aFile) throws FileNotFoundException
	{
		FileOutputStream fos = new FileOutputStream(aFile);
		PrintStream ps = new PrintStream(fos);
		double a = 1;
		ps.print(makeAtomsHeader(a));
		int aUnits = 1;
		
		for(int i = 0; i < bigLattice.length; i++)
		{
			ps.print(bigLattice[i].toStringForAtoms(a, a, a, 0, aUnits, aUnits, aUnits));
		}
		
	}
	public void printLatticeToFile(File aFile, IdealTetrahedron[] lattice) throws FileNotFoundException {
		FileOutputStream fos = new FileOutputStream(aFile);
		PrintStream ps = new PrintStream(fos);
		double a = 1;
		ps.print(makeAtomsHeader(a));
		int aUnits = 1;
		
		for(int i = 0; i < lattice.length; i++)
		{
			ps.print(lattice[i].toStringForAtoms(a, a, a, 0, aUnits, aUnits, aUnits));
		}
		
	}
	public void printLatticeToFile(File aFile, IdealTetrahedron[][] lattice) throws FileNotFoundException {
		FileOutputStream fos = new FileOutputStream(aFile);
		PrintStream ps = new PrintStream(fos);
		ps.println(lattice.length * lattice[0].length);
		ps.println(lattice.length);
		ps.println(lattice[0].length);
		for(int i = 0; i < lattice.length; i++) {
			ps.println(i);
			for(int j = 0; j < lattice[i].length; j++) {
				ps.print(lattice[i][j].toStringForXYZ());
			}
		}
	}
	public void printFooToFile(File aFile) throws FileNotFoundException {
		FileOutputStream fos = new FileOutputStream(aFile);
		PrintStream ps = new PrintStream(fos);
		double a = 1;
		ps.print(makeAtomsHeader(a));
		double scalar = aSmall+10;
		
		// translate the miniLattices
		int y = 0;
		int x = 1;
		for(int i = 0; i < 12; i++) {
			for(int j = 0; j < rotatedEightCubes[i].length; j++) {
				rotatedEightCubes[i][j].translate(new JVector(scalar*x, scalar*y, 0));
			}
			x++;
		}
		
		for(int i = 0; i < 12; i++) {
			for(int j = 0; j < rotatedEightCubes[i].length; j++) {
				ps.print(rotatedEightCubes[i][j].toStringForAtoms(a, a, a, 1, 1, 1, 1));
			}
			               
		}
	}
	public IdealTetrahedron[][] sort(IdealTetrahedron[][] lattice, int whichBox) {
		DoubleLinkedListCBr4 dll = new DoubleLinkedListCBr4(5);
		IdealTetrahedron[][] sorted = new IdealTetrahedron[lattice.length][lattice[0].length];
		sorted = new IdealTetrahedron[lattice.length][lattice[0].length];
		JVector[] sites = JVector.getFCCPositions(whichBox);
		for(int i = 0; i < sites.length; i++) {
			sites[i] = JVector.multiply(sites[i], 4.42);
		}
		double dist;
		for(int i = 0; i < lattice.length; i++) {
			for(int j = 0; j < lattice[0].length; j++) {
				dll.clear();
				for(int k = 0; k < sites.length; k++) {
					dist = JVector.subtract(lattice[i][j].getCenter().getPosition(), sites[k]).length();
					dist = Math.abs(dist);
					dll.insert(new TwoNodeCBr4(lattice[i][j], dist, null, null));
				}
				sorted[i][j] = (IdealTetrahedron) dll.removeHead().getValue().clone();
			}
		}
		return sorted;
	}
	public static void main(String[] args) {
		LennardJonesPotential ljp = new LennardJonesPotential(35, 35, 3.3854, .132);
		
		LennardJonesLookup ljl = new LennardJonesLookup(.01, ljp);
		
		MonoclinicCoordsToCubic m = new MonoclinicCoordsToCubic();
		
		m.aUnits = 2;
		m.aSmall = 2*8.82;
		m.aBig = 55;
		m.aPlastic = 8.82;
		m.aMono = 17.2;
		try {
			//m.readFile(new File("monoclinic crystallographic coords.txt"));
			m.readFile2(new File("monoclinic crystallographic coords.txt"));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}

		double startingSize = 2; // 2x2x2
		double finishSize = 1.5; // 1x1x1
		//double finishSize = 1.5; // 1.5x1.5x1.5
		m.extendCrystal((int) startingSize);
		m.moreMonoCoordsToCartesianCubic();
		m.getCubeFromPosition(uniqueCenters_cartesian[0], (int) startingSize);
		m.makeEightCubes((int) startingSize);
		m.idealizeTetrahedra(1.91);
		m.rotateEightCubesx12();
		m.lookAtDisplacements(m.rotatedEightCubes[0], (int) startingSize);
		m.zeroEightCubes();
		//m.sort(m.rotatedEightCubes, (int) Math.rint(startingSize));
		//m.sort(m.eightCubes, (int) Math.rint(startingSize));
		try {
			m.truncatedEightCubes = m.truncateCubes(m.rotatedEightCubes, (int) (2*startingSize), (int) (2*finishSize));
			m.printLatticeToFile(new File(finishSize + "_truncatedRotatedBoxes.xyz"), m.truncatedEightCubes);
			m.truncatedEightCubes = m.truncateCubes(m.eightCubes, (int) (2*startingSize), (int) (2*finishSize));
			m.printLatticeToFile(new File(finishSize + "_truncatedBoxes.xyz"), m.truncatedEightCubes);
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		
		try {
			m.printLatticeToFile(new File("before_rotation.xyz"), m.eightCubes);
			m.printLatticeToFile(new File("intactBoxes.xyz"), m.rotatedEightCubes);
			m.printLatticeToFile(new File(finishSize + "_truncatedBoxes.xyz"), m.truncatedEightCubes);
			for(int i = 0; i < 8; i++)
				m.printLatticeToFile(new File(i + "test_eight_boxes.inp"), m.eightCubes[i]);
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
	
		try {
			m.printFooToFile(new File("test_rotated_boxes"));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		
		m.rotatePseudoLattice();
		m.truncatePseudoLattices();
		m.getFourSmallCubes();
		m.makeFourBigCubes();
		m.rotateBigCubes();
		m.truncateBigCubes();
		m.makeBigBox(m.aUnits, m.aUnits, m.aUnits);
		/*
		try {
			m.printToFile(new File("1pseudocube.inp"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		*/
		//m.printBigBoxToConsole();
		
		//CBr4PlasticCrystalLattice lattice = new CBr4PlasticCrystalLattice(5, m.aSmall);
		
		//lattice.molecules = m.bigLattice;
		
		//SimulateLattice sl = new SimulateLattice(lattice, ljl);
		
		//sl.buildSurroundings();
		
		try {
			m.printFooToFile(new File("blah.inp"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		//sl.randomWalk(4, 512);
		
		
		
		/*
		try {
			m.makeLatticeFromFinishedBoxFile(new File("monoclinic_cubic_box.inp"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		m.makeBigBox();
		m.printBigBoxToConsole();
		  */
		 
	}
}
