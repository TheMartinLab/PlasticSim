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
public class IdealTetrahedronTest {

	public static void main(String[] args) {
		JAtom center = new JAtom(6, new JVector(0, 0, 0));
		JAtom ligand = new JAtom(35, new JVector(1, 1, 1));
		IdealTetrahedron test = new IdealTetrahedron(center, ligand);
		IdealTetrahedron rotated;
		JVector[] v110s = JVector.v110s;
		JVector trans = new JVector(0, 0, 0);
		for(int i = 0; i < v110s.length; i++) {
			rotated = (IdealTetrahedron) test.clone();
			rotated.align2Fold(v110s[i]);
			for(int j = 2; j < 4; j++) {
				trans = JVector.add(trans, JVector.subtract(rotated.ligands[j].getPosition(), rotated.center.getPosition()));
			}
			trans = JVector.multiply(trans, .25);
			rotated.translate(new JVector(4*i, 0, 0));
			rotated.translate(trans);
			System.out.print(rotated.toStringForXYZ());
			System.out.println(7 + ", " + JVector.add(v110s[i], new JVector(4*i, 0, 0)));	
		}
		
		for(int i = 0; i < v110s.length; i++) {
			System.out.println(v110s[i]);
		}
		
	}
}
