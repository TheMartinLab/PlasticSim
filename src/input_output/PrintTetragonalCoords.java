package input_output;

import defaultPackage.JVector;
import defaultPackage.Lattice;

public class PrintTetragonalCoords {
	public static void main(String[] args) {
		Lattice.printOrientations();
		JVector[][] axes = JVector.axes111U_XYZ;

		JVector[][] axes1 = new JVector[axes.length][axes[0].length];
		axes1[0] = axes[0];
		
		JVector axis = new JVector(0, 0, 1);
		double phi = 90;
		for(int i = 1; i < 4; i++) {
			for(int j = 0; j < axes[0].length; j++) {
				axes1[i][j] = JVector.rotate(axes1[i-1][j], axis, JVector.zero, phi).roundInt();
			}
		}
		
		
		

		for(int i = 0; i < axes.length; i++) {
			for(int j = 0; j < axes[0].length; j++) {
				System.out.print(axes[i][j].toString() + "\t");
			}
			System.out.println();
		}
		axes1 = JVector.axes111_uniqueCrystallographic;
		for(int i = 0; i < axes1.length; i++) {
			for(int j = 0; j < axes1[0].length; j++) {
				System.out.print(axes1[i][j].toString() + "\t");
			}
			System.out.println();
		}
	}
}
