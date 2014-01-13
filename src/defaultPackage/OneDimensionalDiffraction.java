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
import java.text.DecimalFormat;
import java.util.Vector;


public class OneDimensionalDiffraction {
	private ComplexScatteringFactor csf;
	public OneDimensionalDiffraction() {
		csf = new ComplexScatteringFactor();
	}
	public double[] createBins(double qStep, double qMax) {
		double[] bins = new double[(int) Math.rint(qMax / qStep)];
		DecimalFormat twoDForm = new DecimalFormat("#.##");
		for(int i = 0; i < bins.length; i++) {
			bins[i] = Double.valueOf(twoDForm.format(qStep * i));
		}
		twoDForm = null;
		return bins;
	}
	private Vector<Vector<JVector>> initVectorArray(int howMany) {
		Vector<Vector<JVector>> theVectors = new Vector<Vector<JVector>>();
		for(int i = 0; i < howMany; i++) {
			theVectors.add(i, new Vector<JVector>(10));
		}
		return theVectors;
	}
	private Vector<Vector<JPixel>> initPixelArray(int howMany) {
		Vector<Vector<JPixel>> thePixels = new Vector<Vector<JPixel>>();
		for(int i = 0; i < howMany; i++) {
			thePixels.add(i, new Vector<JPixel>(10));
		}
		return thePixels;
	}
	public Vector<Vector<JVector>> getScatteringVectors(double qStep, double qMax) {
		int max = (int) Math.rint(qMax/qStep);
		int min = -max;
		Vector<Vector<JVector>> theVectors = initVectorArray(max+1);
		DecimalFormat twoDForm = new DecimalFormat("#.###");
		JVector cur = new JVector();
		double len;
		int binIdx;
		for(int h = min; h <= max; h++) {
			cur.i = Double.valueOf(twoDForm.format(h*qStep));
			for(int k = min; k <= max; k++) {
				cur.j = Double.valueOf(twoDForm.format(k*qStep));
				for(int l = min; l <= max; l++) {
					cur.k = Double.valueOf(twoDForm.format(l*qStep));
					len = cur.length();
					if(len < qMax && len != 0) {
						binIdx = (int) (len / qStep);
						theVectors.get(binIdx).add((JVector) cur.clone());
					}
				}
			}
		}
		return theVectors;
	}
	public Vector<Vector<JPixel>> createAllPixels(int[] elemTypes, double qStep, double qMax, double wavelength) throws Exception {
		Vector<Vector<JVector>> theVectors = getScatteringVectors(qStep, qMax);
		Vector<Vector<JPixel>> thePixels = new Vector<Vector<JPixel>>();
		for(int i = 0; i < theVectors.size(); i++) {
			thePixels.add(createPixels(theVectors.get(i), elemTypes, wavelength));
		}
		return thePixels;
	}
	public Vector<JPixel> createPixels(Vector<JVector> theVectors, int[] elemTypes, double wavelength) throws Exception {
		Vector<JPixel> thePixels = new Vector<JPixel>();
		JVector curVec;
		JPixel curPix;
		JComplex sf;
		for(int i = 0; i < theVectors.size(); i++) {
			curVec = theVectors.get(i);
			curPix = new JPixel(elemTypes.length, (JVector) curVec.clone(), elemTypes);
			for(int elemTypeIdx = 0; elemTypeIdx < elemTypes.length; elemTypeIdx++) {
				sf = csf.generateF0(curVec, elemTypes[elemTypeIdx], csf.getComplexF(1, wavelength, elemTypes[elemTypeIdx]));
				curPix.setSf(elemTypeIdx, sf);
			}
			thePixels.add(curPix);
		}
		return thePixels;
	}

}
