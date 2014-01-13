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
import java.util.Random;
import java.io.*;
import java.util.*;

public class IdealTetrahedron extends Molecule implements Serializable 
{
	private static final long serialVersionUID = 5257594070437606835L;

	JAtom center;
	
	JAtom[] ligands;
	
	int orientation;
	
	private static int numTotal;
	
	private final int ID;
	
	//protected IdealTetrahedron[] surrounding;

	/**
	 * This constructor initializes an ideal tetrahedron with ligands in the four 111, -1 -1 1, -1 1 -1, 1 -1 -1 
	 * directions and angles between the atoms = 109.47
	 * @param center	The central atom of the tetrahedron
	 * @param ligands	The four ligands of the tetrahedron
	 */
	public IdealTetrahedron(JAtom center, JAtom aLigand)
	{
		numTotal++;
		
		ID = numTotal;
		this.center = center;
		
		ligands = new JAtom[4];
		
		ligands[0] = aLigand;

		//surrounding = new IdealTetrahedron[12];
		
		generateLigands();
	}
	public JAtom[] getAtoms() {
		JAtom[] temp = new JAtom[5];
		temp[0] = center;
		temp[1] = ligands[0];
		temp[2] = ligands[1];
		temp[3] = ligands[2];
		temp[4] = ligands[3];
		return temp;
	}
	public IdealTetrahedron(JAtom center, JAtom[] ligands)
	{
		numTotal++;
		
		ID = numTotal;
		
		this.center = center;
		
		this.ligands = ligands;
	}
	public IdealTetrahedron(JAtom[] atoms)
	{
		numTotal++;
		
		ID = numTotal;
		
		this.center = (JAtom) atoms[0].clone();
		this.ligands = new JAtom[4];
		for(int i = 0; i < 4; i++) {
			this.ligands[i] = (JAtom) atoms[i+1].clone();
		}
		
	}
	public void generateLigands()
	{
		/* Define the origin */
		JVector origin = (JVector) center.getPosition().clone();
		
		/* Define the axes */
		JVector x = new JVector(1, 0, 0);
		JVector y = new JVector(0, 1, 0);
		JVector z = new JVector(0, 0, 1);
		JVector[] axis = {x, y, z};

		/* Define the angle */
		double phi = 180;
		
		/* Convert an JAtom into a Quaternion */
		Quaternion temp = null;
		
		/* Rotate the JAtom into the other three positions */
		for(int i = 1; i < 4; i++)
		{
			temp = new Quaternion(ligands[0].getPosition());
			
			temp = Quaternion.rotate(temp, axis[i-1], origin, phi);
			
			ligands[i] = new JAtom(ligands[0].Z, temp.position);
		}
	}

	@Override
	public void rotate(JVector axis, JVector origin, double phi) 
	{
		Quaternion temp = new Quaternion(center.getPosition());
		
		temp = Quaternion.rotate(temp, axis, origin, phi);
		
		center.setPosition(temp.getIJK());
		
		for(int i = 0; i < 4; i++)
		{
			temp = new Quaternion(ligands[i].getPosition());
			
			temp = Quaternion.rotate(temp, axis, origin, phi);
			
			ligands[i].setPosition(temp.getIJK());
		}
		
	}

	@Override
	public void moveTo(JVector newPos) {
		JVector translationVec = JVector.subtract(newPos, center.getPosition());
		
		center.setPosition(translationVec);
		
		for(int i = 0; i < ligands.length; i++) {
			ligands[i].translate(translationVec);
		}
	}
	@Override
	public void translate(JVector amount)
	{
		center.translate(amount);
		
		for(int i = 0; i < 4; i++)
		{
			ligands[i].translate(amount);
		}
	}
	public void align2Fold(JVector toAlign) {
		JVector curPos = JVector.add(JVector.subtract(ligands[0].getPosition(), center.getPosition()), JVector.subtract(ligands[1].getPosition(), center.getPosition()));
		double phi = JVector.angle(toAlign, curPos);
		if(phi == 0)
			return;
		while(phi == 180) {
			perturb(.01);
			curPos = JVector.add(JVector.subtract(ligands[0].getPosition(), center.getPosition()), JVector.subtract(ligands[1].getPosition(), center.getPosition()));
			phi = JVector.angle(toAlign, curPos);
		}
		rotate(JVector.cross(toAlign, curPos), center.getPosition(), -phi);
	}
	public void alignBr(JVector toAlign) {
		double phi = JVector.angle(toAlign, JVector.subtract(ligands[0].getPosition(), center.getPosition()));
		if(phi == 0)
			return;
		while(phi == 180) {
			perturb(.01);
			phi = JVector.angle(toAlign, JVector.subtract(ligands[0].getPosition(), center.getPosition()));
		}
		rotate(JVector.cross(toAlign, JVector.subtract(ligands[0].getPosition(), center.getPosition())), center.getPosition(), -phi);
	}
	public void perturb(double phiToPerturb) {
		JVector axis = new JVector(Math.random(), Math.random(), Math.random());
		rotate(axis, center.getPosition(), phiToPerturb);
	}
	
	public void align2ndBr(JVector firstAlignment, int whichToAlignAlong) {
		String idx = "";
		double phi;
		for(int i = 0; i < JVector.v110s.length; i++) {
			phi = JVector.angle(firstAlignment, JVector.v110s[i]);
			if(phi >= 90 && phi <= 120)
				idx += i + "_";
		}
		String[] strChoices = idx.split("_");
		int choice = Integer.valueOf(strChoices[whichToAlignAlong]);
		// rotate the chosen vector so that it is 109.47° away from the initial alignment
		JVector target = (JVector) JVector.v110s[choice].clone();
		JVector axis = JVector.cross(target, firstAlignment);
		phi = JVector.angle(target, firstAlignment);
		target = JVector.rotate(target, axis, new JVector(0, 0, 0), phi - 109.47);
		phi = JVector.angle(target, firstAlignment);
		// choose a bromine to align
		phi = JVector.angle(target, JVector.subtract(ligands[1].getPosition(), center.getPosition()));
		int whileIdx = 0;
		while(Math.abs(phi) > .01) {
			rotate(firstAlignment, center.getPosition(), -phi);
			phi = JVector.angle(target, JVector.subtract(ligands[1].getPosition(), center.getPosition()));
			whileIdx++;
			if(whileIdx > 100) {
				return;
			}
		}
	}
	/**
	 * @return  the vector position of the center(0) and the four ligands(1-4)
	 */
	public JVector[] getVectors()
	{
		JVector atoms[] = new JVector[5];
		
		atoms[0] = center.getPosition();
		
		for(int i = 1; i < 5; i++)
		{
			atoms[i] = ligands[i-1].getPosition();
		}
		
		return atoms;
	}
	
	/**
	 * Multiply all the positions by an amount, a
	 * @param a The amount to multiply by
	 */
	public void scaleByAmount(double a)
	{
		center.setPosition(JVector.multiply(center.getPosition(), a));
		for(int i = 0; i < ligands.length; i++)
		{
			ligands[i].setPosition(JVector.multiply(ligands[i].getPosition(), a));
		}
	}

	public JAtom[] getLigands()
	{
		JAtom[] cloned = new JAtom[ligands.length];
		
		for(int i = 0; i < ligands.length; i++)
		{
			cloned[i] = (JAtom) ligands[i].clone();
		}
		
		return cloned;
	}
	
	public void setLigands(JAtom[] newLigands) { ligands = newLigands; }
	
	//public void setSurrounding(IdealTetrahedron[] surrounding) { this.surrounding = surrounding; }
	//public IdealTetrahedron[] getSurrounding() { return surrounding; }
	public String toString()
	{
		String ligands = "";
		
		for(int i = 0; i < 4; i++)
		{
			ligands += this.ligands[i].atomID + "," + this.ligands[i] + "\n";
		}
		return center.atomID + "," + center + "\n" + ligands;
	}

	@Override
	public String toStringForAtoms(double a, double b, double c, int number, int aUnits, int bUnits, int cUnits) 
	{
		String s = center.toStringForAtoms(a, a, a, center.atomID, aUnits, bUnits, cUnits) + "\n";
		
		for(int i = 0; i < 4; i++)
		{
			int number2 = number + 1;
			if(i != 0)
				number2 = 4*number + i;
			s+= ligands[i].toStringForAtoms(a, a, a, ligands[i].atomID, aUnits, bUnits, cUnits) + "\n";
		}
		
		return s;
	}
	
	public String toStringForXYZ() 
	{
		String s = center.toStringForXYZ();
		
		for(int i = 0; i < 4; i++)
		{
			s+= ligands[i].toStringForXYZ();
		}
		
		return s;
	}
	
	public Object clone()
	{
		JAtom tempCenter = (JAtom) center.clone();

		JAtom[] tempLigands = new JAtom[4];
		
		for(int i = 0; i < 4; i++)
		{
			tempLigands[i] = (JAtom) ligands[i].clone();
		}
		return new IdealTetrahedron(tempCenter, tempLigands);
	}
	
	public JAtom getCenter()
	{
		return center;
	}
	
	public void setCenter(JVector toSet)
	{
		center.getPosition().i = toSet.i;
		
		center.getPosition().j = toSet.j;
		
		center.getPosition().k = toSet.k;
	}
	
	public int getID() { return ID; }

	public Molecule[] getSurrounding() {
		return surrounding;
	}
	
	public void setSurrounding(IdealTetrahedron[] surroundingCBr4) {
		surrounding = surroundingCBr4;
	}
}
