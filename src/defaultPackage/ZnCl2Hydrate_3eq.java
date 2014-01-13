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
public class ZnCl2Hydrate_3eq {

	private double a = 6.55;
	private double b = 6.6;
	private double c = 15.12;
	private static boolean inversionSwitch = false;
	private final static double ZnCl_len = 2.24;
	private final static double ZnO_len = 2.04;
	public static Molecule getHexaAquoZinc() {
		JAtom zn = new JAtom(30, new JVector(0, 0, 0));
		IdealOctahedron hexaAquoZinc = new IdealOctahedron(zn, 8, ZnO_len);
		return hexaAquoZinc;
	}
	public static Molecule getTetrachloroZincate() {
		if(inversionSwitch) {
			inversionSwitch = !inversionSwitch;
			return getTetrachloroZincateOne();
		}
		
		else {
			inversionSwitch = !inversionSwitch;
			return getTetrachloroZincateTwo();
		}
	}
	/**
	 * Inversion of Two
	 * @return
	 */
	public static Molecule getTetrachloroZincateOne() {
		JAtom zn = new JAtom(30, new JVector(0, 0, 0));
		double bondLen = ZnCl_len/Math.sqrt(3);
		JAtom cl = new JAtom(17, new JVector(bondLen, bondLen, bondLen));
		IdealTetrahedron tetraChloroZincate = new IdealTetrahedron(zn, cl);
		return tetraChloroZincate;
	}
	/**
	 * Inversion of One
	 * @return
	 */
	public static Molecule getTetrachloroZincateTwo() {
		IdealTetrahedron tetraChloroZincate = (IdealTetrahedron) getTetrachloroZincate();
		JAtom[] ligands = tetraChloroZincate.ligands;
		JVector newPos;
		for(int i = 0; i < ligands.length; i++) {
			newPos = JVector.multiply(ligands[i].getPosition(), -1);
			ligands[i].setPosition(newPos);
		}
		return tetraChloroZincate;
	}
	public static Molecule[] getBasisMolecules() {
		Molecule[] basis = new Molecule[3];
		
		basis[0] = getHexaAquoZinc();
		basis[0].rotate(JVector.z, basis[0].getCenter().getPosition(), 23);
		
		basis[1] = getTetrachloroZincate();
		((IdealTetrahedron) basis[1]).rotate(JVector.x, basis[1].getCenter().getPosition(), -45);
		((IdealTetrahedron) basis[1]).rotate(JVector.z, basis[1].getCenter().getPosition(), -35.25);
		
		basis[2] = (Molecule) basis[1].clone();
		((IdealTetrahedron) basis[2]).rotate(JVector.x, basis[1].getCenter().getPosition(), 180);
		((IdealTetrahedron) basis[2]).rotate(JVector.y, basis[1].getCenter().getPosition(), 180);
		
		return basis;
	}
	public static Molecule[] getUnitCellMolecules() {
		Molecule[] basis = new Molecule[4];
		
		basis[0] = getHexaAquoZinc();
		basis[0].rotate(JVector.z, basis[0].getCenter().getPosition(), 23);
		basis[1] = (Molecule) basis[0].clone();
		
		basis[2] = getTetrachloroZincate();
		((IdealTetrahedron) basis[2]).rotate(JVector.x, basis[1].getCenter().getPosition(), -45);
		((IdealTetrahedron) basis[2]).rotate(JVector.z, basis[1].getCenter().getPosition(), -35.25);
		
		basis[3] = (Molecule) basis[2].clone();
		((IdealTetrahedron) basis[3]).rotate(JVector.x, basis[1].getCenter().getPosition(), 180);
		((IdealTetrahedron) basis[3]).rotate(JVector.y, basis[1].getCenter().getPosition(), 180);
		
		basis[0].translate(new JVector(3.275, 3.3, 0));
		basis[1].translate(new JVector(3.275, 3.3, 7.56));		
		basis[2].translate(new JVector(1.070925, 5.841, 3.78));
		basis[3].translate(new JVector(5.479075, .759, 11.34));
		return basis;
	}
	public static void main(String[] args) {
		
		ZnCl2Hydrate_3eq test = new ZnCl2Hydrate_3eq();
		Molecule[] basis = getUnitCellMolecules();
		
		System.out.println(basis[0]);
		System.out.println(basis[1]);
		System.out.println(basis[2]);
		System.out.println(basis[3]);
	}
}
