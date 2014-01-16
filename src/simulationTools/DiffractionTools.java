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

import java.util.HashMap;

import defaultPackage.ComplexScatteringFactor;
import defaultPackage.JComplex;
import defaultPackage.JVector;
import defaultPackage.JPixel;

public class DiffractionTools {
	/**
	 * 
	 * @param vFamily in order: XYZ
	 * @param qxMax
	 * @param qyMax
	 * @param qStep
	 * @param elemTypes
	 * @param wavelength
	 * @return
	 */
	public static JPixel[] initPixels(JVector[][] vFamily, double qxMax, double qyMax, 
			double qStep, int[] elemTypes, double wavelength) {
		int numPixels = 0;
		int pixIdx;
		JVector Q, vx, vy;
		JComplex sf;

		for(int i = 0; i < vFamily.length; i++) {
			for(double qy = -qyMax; qy <= qyMax; qy += qStep) {
				for(double qx = -qxMax; qx <= qxMax; qx += qStep) {
					numPixels++;
				}
			}
		}
		JPixel[] p = new JPixel[numPixels];
		System.out.println("numPixels: " + numPixels);
		ComplexScatteringFactor csf = new ComplexScatteringFactor();
		pixIdx = 0;
		for(pixIdx = 0; pixIdx < p.length; pixIdx++) { p[pixIdx] = new JPixel(elemTypes.length); }
		pixIdx = 0;
		/* init new pixels */
		for(int elemTypeIdx = 0; elemTypeIdx < elemTypes.length; elemTypeIdx++) {
			pixIdx = 0;
			for(int i = 0; i < vFamily.length; i++) {
				vx = vFamily[i][0];
				vy = vFamily[i][1];
				HashMap<Double, JComplex> sfMap = csf.buildHashMap(Math.max(qxMax, qyMax), qStep, wavelength, elemTypes[elemTypeIdx], vx, vy);
				for(double qy = -qyMax; qy <= qyMax; qy += qStep) {
					for(double qx = -qxMax; qx <= qxMax; qx += qStep) {
						Q = JVector.add(JVector.multiply(vx, qx), JVector.multiply(vy, qy));
						sf = sfMap.get(Q.length());
						
						if(sf == null) {
							try { 
								sf = csf.generateF0(Q, elemTypes[elemTypeIdx], csf.getComplexF(1, wavelength, elemTypes[elemTypeIdx])); 
								//System.out.println("elemTypes[elemTypeIdx]: " + elemTypes[elemTypeIdx]);
							} catch(Exception e) { e.printStackTrace(); }
						}
						//System.out.println(qx + "," + qy+ "," + sf);
						if(elemTypeIdx == 0) { 
							p[pixIdx].setQ((JVector) Q.clone());
						}
						p[pixIdx].setAtomTypes(elemTypes);
						p[pixIdx].setSf(elemTypeIdx, (JComplex) sf.clone());
						p[pixIdx].setI(0);
						pixIdx++;
					}
				}
			}
		}

		return p;
	}
	/**
	 * 
	 * @param vFamily in order: ZXY
	 * @param qxMax
	 * @param qyMax
	 * @param qStep
	 * @param elemTypes
	 * @param wavelength
	 * @return
	 */
	public static JPixel[] initPixels_old(JVector[][] vFamily, double qxMax, double qyMax, 
			double qStep, int[] elemTypes, double wavelength, double zShift) {
		int numPixels = 0;
		int pixIdx;
		JVector Q, vx, vy, vz;
		JComplex sf;

		for(int i = 0; i < vFamily.length; i++) {
			for(double qy = -qyMax; qy <= qyMax; qy += qStep) {
				for(double qx = -qxMax; qx <= qxMax; qx += qStep) {
					numPixels++;
				}
			}
		}
		JPixel[] p = new JPixel[numPixels];
		System.out.println("numPixels: " + numPixels);
		ComplexScatteringFactor csf = new ComplexScatteringFactor();
		pixIdx = 0;
		for(pixIdx = 0; pixIdx < p.length; pixIdx++) { p[pixIdx] = new JPixel(elemTypes.length); }
		pixIdx = 0;
		/* init new pixels */
		for(int elemTypeIdx = 0; elemTypeIdx < elemTypes.length; elemTypeIdx++) {
			pixIdx = 0;
			for(int i = 0; i < vFamily.length; i++) {
				vx = vFamily[i][1];
				vy = vFamily[i][2];
				vz = vFamily[i][0];
				HashMap<Double, JComplex> sfMap = csf.buildHashMap(Math.max(qxMax, qyMax), qStep, wavelength, elemTypes[elemTypeIdx], vx, vy);
				for(double qy = -qyMax; qy <= qyMax; qy += qStep) {
					for(double qx = -qxMax; qx <= qxMax; qx += qStep) {
						Q = JVector.add(JVector.multiply(vx, qx), JVector.multiply(vy, qy));
						Q = JVector.add(Q, JVector.multiply(vz, zShift));
						sf = sfMap.get(Q.length());
						
						if(sf == null) {
							try { 
								sf = csf.generateF0(Q, elemTypes[elemTypeIdx], csf.getComplexF(1, wavelength, elemTypes[elemTypeIdx])); 
								//System.out.println("elemTypes[elemTypeIdx]: " + elemTypes[elemTypeIdx]);
							} catch(Exception e) { e.printStackTrace(); }
						}
						//System.out.println(qx + "," + qy+ "," + sf);
						if(elemTypeIdx == 0) { 
							p[pixIdx].setQ((JVector) Q.clone());
						}
						p[pixIdx].setSf(elemTypeIdx, (JComplex) sf.clone());
						p[pixIdx].setI(0);
						pixIdx++;
					}
				}
			}
		}

		return p;
	}
}
