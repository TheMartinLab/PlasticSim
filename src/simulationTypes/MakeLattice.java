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
package simulationTypes;

import java.util.Observable;
import java.util.Observer;

import defaultPackage.Lattice;

public class MakeLattice extends Observable {

	int numUnits = 10;
	double a = 8.82;
	private Lattice l;
	
	public MakeLattice() {
		
	}
	
	
}
