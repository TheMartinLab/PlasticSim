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
import io.MyObjectInputStream;
import io.MyPrintStream;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Vector;

import chemistry.JAtomTools;
import simulationTypes.Lattices;
import defaultPackage.JAtom;
import defaultPackage.JVector;
import defaultPackage.Molecule;


public class LatticeTools {
	public static void printLatticeXYZ(JAtom[] lattice, PrintStream ps) {
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
	}
	public static void printLatticeXYZ(Molecule[][][] lattice, PrintStream ps, double a) {
		printLatticeXYZ(molToFractionalAtoms(lattice, a), ps);
	}
	public static JAtom[] molToFractionalAtoms(Molecule[][][] lattice, double a) {
		Vector<JAtom> atoms = new Vector<JAtom>();
		JAtom curAtom;
		JAtom[] mol;
		JVector curVector;
		
		for(int i = 0; i < lattice.length; i++){
			for(int j = 0; j < lattice[i].length; j++){
				for(int k = 0; k < lattice[i][j].length; k++){
					if(lattice[i][j][k] != null) {
						mol = lattice[i][j][k].getAtoms();
						for(int atomIdx = 0; atomIdx < mol.length; atomIdx++) {
							curAtom = (JAtom) mol[atomIdx].clone();
							curVector = curAtom.getPosition();
							curVector.multiply(1./a);
							atoms.add(curAtom);
						}
					}
				}
			}
		}
		
		return atoms.toArray(new JAtom[atoms.size()]);
	}
	public static Molecule[] _3dTo1d(Molecule[][][] lattice) {
		Vector<Molecule> molecules = new Vector<Molecule>();
		for(int a = 0; a < lattice.length; a++) {
			for(int b = 0; b < lattice[a].length; b++) {
				for(int c = 0; c < lattice[a][b].length; c++) {
					if(lattice[a][b][c] != null)
						molecules.add(lattice[a][b][c]);
				}
			}
		}
		
		return molecules.toArray(new Molecule[molecules.size()]);
	}
	public static JAtom[] moleculesToAtoms(Molecule[][][] lattice) {
		Vector<JAtom> atoms = new Vector<JAtom>();
		JAtom[] mol;
		
		for(int i = 0; i < lattice.length; i++){
			for(int j = 0; j < lattice[i].length; j++){
				for(int k = 0; k < lattice[i][j].length; k++){
					if(lattice[i][j][k] != null) {
						mol = lattice[i][j][k].getAtoms();
						for(int atomIdx = 0; atomIdx < mol.length; atomIdx++) {
							atoms.add(mol[atomIdx]);
						}
					}
				}
			}
		}
		
		return atoms.toArray(new JAtom[atoms.size()]);
	}
	public static synchronized Lattices readInLattice(File objFile) throws ClassNotFoundException, IOException {
		MyObjectInputStream mois = new MyObjectInputStream(objFile);
		
		Lattices init = (Lattices) mois.readObject();
		
		mois.close();
		
		return init;
	}
	
	private static int zToIdx(int Z) {
		if(Z == 6)
			return 1;
		return 2;
	}
	public static void LatticeToLAMMPS_Atomic(File outputFile, File inputObjFile) {
		Lattices init = null;
		try {
			init = LatticeTools.readInLattice(inputObjFile);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		MyPrintStream mps = new MyPrintStream(outputFile, false);
		
		Molecule[] lattice = _3dTo1d(init.getLattice());
		
		double 	xMin = Double.MAX_VALUE, 
				xMax = Double.MIN_VALUE, 
				yMin = Double.MAX_VALUE, 
				yMax = Double.MIN_VALUE, 
				zMin = Double.MAX_VALUE, 
				zMax = Double.MIN_VALUE,
				x, y, z;

		Vector<String> atomLines = new Vector<String>();
		int atomIdx = 1;
		for(Molecule mol : lattice) {
			for(JAtom atom : mol.getAtoms()) {
				JVector atomPos = atom.getPosition();
				x = atomPos.getI();
				y = atomPos.getJ();
				z = atomPos.getK();
				
				if(xMin > x) { xMin = x; }
				if(yMin > y) { yMin = y; }
				if(zMin > z) { zMin = z; }
				
				if(xMax < x) { xMax = x; }
				if(yMax < y) { yMax = y; }
				if(zMax < z) { zMax = z; }
				
				String atomLine = atomIdx + " " + 
						zToIdx(atom.getZ()) + " " + 
						x + " " +
						y + " " +
						z;
				atomLines.add(atomLine);
				
				atomIdx++;
			}
		}
		
		/* ********************** */
		/* WRITING THE INPUT FILE */
		/* ********************** */
		mps.println("atom_style atomic");
		
		mps.println();
		
		mps.println(lattice.length * 5 + " atoms");
		
		mps.println();
		
		mps.println("2 atom types");
		
		mps.println();
		
		mps.println(xMin + " " + xMax + " xlo xhi");
		mps.println(yMin + " " + yMax + " ylo yhi");
		mps.println(zMin + " " + zMax + " zlo zhi");
		
		mps.println();
		
		mps.println("Masses");
		
		mps.println();
		
		mps.println("1 " + JAtomTools.getMass(6));
		mps.println("2 " + JAtomTools.getMass(35));
		
		
		mps.println("\nAtoms\n");
		int lineIdx = 1;
		String lines = "";
		int printEvery = 1000;
		for(String line : atomLines) {
			if(lineIdx % printEvery == 0) {
				lineIdx = 0;
				mps.print(lines);
				lines = "";
			}
			lines += line + "\n";
			lineIdx++;
		}
		mps.print(lines);
		mps.flush();
		mps.close();
		
	}
	public static void LatticeToLAMMPS_Molecular(File outputFile, File inputObjFile) {
		Lattices init = null;
		try {
			init = LatticeTools.readInLattice(inputObjFile);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		MyPrintStream mps = new MyPrintStream(outputFile, false);
		
		Molecule[] lattice = _3dTo1d(init.getLattice());
		
		double 	xMin = Double.MAX_VALUE, 
				xMax = Double.MIN_VALUE, 
				yMin = Double.MAX_VALUE, 
				yMax = Double.MIN_VALUE, 
				zMin = Double.MAX_VALUE, 
				zMax = Double.MIN_VALUE,
				x, y, z;

		Vector<String> atomLines = new Vector<String>();
		Vector<String> bondLines = new Vector<String>();
		int moleculeIdx = 0;
		int bondIdx = 1;
		int atomIdx = 0;
		int bondType;
		/* SET ALL ATOM INDICES */
		for(Molecule mol : lattice)
			for(JAtom atom : mol.getAtoms())
				atom.setAtomID(++atomIdx);

		/* LOOP THROUGH ALL MOLECULES */
		for(Molecule mol : lattice) {
			JAtom[] atoms = mol.getAtoms();
			JAtom atom1, atom2;
			int Z1;
			++moleculeIdx;

			/* LOOP THROUGH ATOMS ON A MOLECULE */
			for(int i = 0; i < atoms.length; i++) {
				atom1 = atoms[i];
				Z1 = atom1.getZ();
				JVector atomPos = atom1.getPosition();
				x = atomPos.getI();
				y = atomPos.getJ();
				z = atomPos.getK();
				
				if(xMin > x) { xMin = x; }
				if(yMin > y) { yMin = y; }
				if(zMin > z) { zMin = z; }
				
				if(xMax < x) { xMax = x; }
				if(yMax < y) { yMax = y; }
				if(zMax < z) { zMax = z; }

				atomLines.add(atom1.getAtomID() + " " + moleculeIdx + " " + zToIdx(Z1) + " " + x + " " + y + " " + z);
				
				/* SET UP INTRAMOLECULAR BONDS */
				for(int j = i; j < atoms.length; j++) {
					atom2 = atoms[j];
					if(atom1 == atom2)
						continue;

					bondType = 2;
					if(Z1 == 6)
						bondType = 1;
					
					bondLines.add(bondIdx++ + " " + bondType + " " + atom1.getAtomID() + " " + atom2.getAtomID());
				}
			}
		}
		
		/* ********************** */
		/* WRITING THE INPUT FILE */
		/* ********************** */
		mps.println("atom_style molecular");
		
		mps.println();
		
		mps.println(lattice.length * 5 + " atoms");
		mps.println(lattice.length * 10 + " bonds");
		
		mps.println();
		
		mps.println("2 atom types");
		mps.println("2 bond types");
		
		mps.println();
		
		mps.println(xMin + " " + xMax + " xlo xhi");
		mps.println(yMin + " " + yMax + " ylo yhi");
		mps.println(zMin + " " + zMax + " zlo zhi");
		
		mps.println();
		
		mps.println("Masses");
		
		mps.println();
		
		mps.println("1 " + JAtomTools.getMass(6));
		mps.println("2 " + JAtomTools.getMass(35));
		
		
		mps.println("\nAtoms\n");
		mps.printAll(atomLines);
		mps.flush();
		
		mps.println("\nBonds\n");
		mps.printAll(bondLines);
		mps.flush();
		mps.close();
	}
}
