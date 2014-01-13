package defaultPackage;

import io.MyPrintStream;

import java.io.File;
import java.util.Vector;

public class MakeLatticeAndRotateAround110 {

	public static void main(String[] args) {
		Vector<IdealTetrahedron> orientations = new Vector<IdealTetrahedron>();
		for(int i = 0; i < 6; i++) {
			orientations.add(Lattice.getOrientation(i));
		}
		
		Vector<JAtom> brominePositions = new Vector<JAtom>();
		
		for(IdealTetrahedron tetra : orientations)
			for(JAtom atom : tetra.getAtoms())
				if(atom.Z == 35)
					brominePositions.add(atom);
		
		brominePositions.add(0, new JAtom(6, new JVector(0, 0, 0)));
		Molecule basisMol = new GenericMolecule(brominePositions);
		
		JVector translation = new JVector();
		
		Vector<Molecule> lattice = new Vector<Molecule>();
		int numUnits = 2;
		double a = 8.82;
		for(int i = 0; i <= numUnits*2; i++) {
			translation.i = i * a / 2.;
			for(int j = 0; j <= numUnits*2; j++) {
				translation.j = j * a / 2.;
				for(int k = 0; k <= numUnits*2; k++) {
					translation.k = k * a / 2.;
					if((i + j + k ) % 2 == 0) {
						Molecule mol = (GenericMolecule) basisMol.clone();
						mol.translate(translation);
						lattice.add(mol);
					}
				}
			}
		}
		
		JVector origin = new JVector(numUnits * a/2., numUnits * a/2., numUnits * a/2.);
		origin.multiply(-1);
		
		for(Molecule mol : lattice)
			mol.translate(origin);
		
		origin = JVector.zero;
		JVector axis = new JVector(1, 1, 0);
		
		double phi = 90;
		
		// perform the twinning operation, a 90 degree rotation about the 110
		int numAtoms = 0; 
		for(Molecule mol : lattice) {
			numAtoms += mol.getAtoms().length;
			mol.rotate(axis, origin, phi);
		}
		
		// now align the 111 with the 001
		JVector v111 = new JVector(1, 1, 1);
		
		double phi1 = 45;
		JVector axis1 = new JVector(0, 0, 1);
		
		JVector v111_2 = JVector.rotate(v111, axis1, JVector.zero, phi1);
		
		double phi2 = JVector.angle(v111, JVector.z);
		JVector axis2 = JVector.x;
		

		for(Molecule mol : lattice) {
			mol.rotate(axis1, origin, phi1);
			mol.rotate(axis2, origin, phi2);
		}
		
		// remove all molecules whose center of mass is +/- 2 angstroms from the xy0 plane
		int numAtomsInSlice = 0;
		Vector<Molecule> slice = new Vector<Molecule>();
		for(Molecule mol : lattice) {
			if(Math.abs(mol.getAtoms()[0].getPosition().k) < a/3.) {
				slice.add(mol);
				numAtomsInSlice += mol.getAtoms().length;
			}
		}
		
		MyPrintStream mps = new MyPrintStream(new File("rotateLatticeTest.xyz"));
		
		mps.println(numAtomsInSlice + "\n");
		
		for(Molecule mol : slice)
			mps.print(mol.toStringForXYZ());
		
		mps.close();
	}
}
