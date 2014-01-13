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

public class IdealOctahedron extends Molecule {

	private JAtom[] ligands;
	private JAtom center;
	
	public IdealOctahedron(JAtom center, int ligandType, double ligandDistance) {
		this.center = center;
		setupLigands(ligandType, ligandDistance);
	}
	public IdealOctahedron(JAtom center, JAtom[] ligands) {
		this.center = center;
		this.ligands = ligands;
	}
	
	private void setupLigands(int ligandType, double ligandDistance) {
		ligands = new JAtom[6];
		JVector pos;
		JAtom atom;
		for(int i = 0; i < JVector.v100s.length; i++) {
			pos = JVector.multiply(JVector.v100s[i], ligandDistance);
			atom = new JAtom(ligandType, pos);
			ligands[i] = atom;
		}
	}
	@Override
	public void rotate(JVector axis, JVector origin, double phi) {
		Quaternion temp = new Quaternion(center.getPosition());
		temp = Quaternion.rotate(temp, axis, origin, phi);
		center.setPosition(temp.getIJK());
		for(int i = 0; i < ligands.length; i++) {
			temp = new Quaternion(ligands[i].getPosition());
			temp = Quaternion.rotate(temp, axis, origin, phi);
			ligands[i].setPosition(temp.getIJK());
		}
	}

	@Override
	public void translate(JVector amount) {
		center.translate(amount);
		for(int i = 0; i < ligands.length; i++) {
			ligands[i].translate(amount);
		}
	}

	public String toString() {
		String ligands = "";
		
		for(int i = 0; i < this.ligands.length; i++) {
			ligands += this.ligands[i].atomID + "," + this.ligands[i] + "\n";
		}
		return center.atomID + "," + center + "\n" + ligands;
	}

	@Override
	public String toStringForAtoms(double a, double b, double c, int number, int aUnits, int bUnits, int cUnits) 	{
		String s = center.toStringForAtoms(a, a, a, center.atomID, aUnits, bUnits, cUnits) + "\n";
		
		for(int i = 0; i < ligands.length; i++) {
			int number2 = number + 1;
			if(i != 0)
				number2 = 4*number + i;
			s+= ligands[i].toStringForAtoms(a, a, a, ligands[i].atomID, aUnits, bUnits, cUnits) + "\n";
		}
		
		return s;
	}
	
	public String toStringForXYZ()  {
		String s = center.toStringForXYZ();
		
		for(int i = 0; i < ligands.length; i++) {
			s+= ligands[i].toStringForXYZ();
		}
		
		return s;
	}
	@Override
	public JAtom[] getLigands() {
		return ligands;
	}
	@Override
	public JAtom getCenter() {
		return center;
	}
	@Override
	public void setLigands(JAtom[] ligands) {
		this.ligands = ligands;
	}
	public JVector[] getVectors() {
		JVector atoms[] = new JVector[1+ligands.length];
		
		atoms[0] = center.getPosition();
		
		for(int i = 1; i < atoms.length; i++) {
			atoms[i] = ligands[i-1].getPosition();
		}
		
		return atoms;
	}
	public Object clone() {
		JAtom newCenter = (JAtom) center.clone();
		JAtom[] newLigands = new JAtom[ligands.length];
		for(int i = 0; i < ligands.length; i++) {
			newLigands[i] = (JAtom) ligands[i].clone();
		}
		return new IdealOctahedron(newCenter, newLigands);
	}
	
	public JAtom[] getAtoms() {
		JAtom[] temp = new JAtom[1+ligands.length];
		temp[0] = center;
		temp[1] = ligands[0];
		temp[2] = ligands[1];
		temp[3] = ligands[2];
		temp[4] = ligands[3];
		temp[5] = ligands[4];
		temp[6] = ligands[5];
		return temp;
	}
	
	@Override
	public void moveTo(JVector newPos) {
		JVector translationVec = JVector.subtract(newPos, center.getPosition());
		
		center.setPosition(translationVec);
		
		for(int i = 0; i < ligands.length; i++) {
			ligands[i].translate(translationVec);
		}
	}
}
