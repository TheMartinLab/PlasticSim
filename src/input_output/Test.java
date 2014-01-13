package input_output;

import defaultPackage.JVector;

public class Test {

	public static void main(String[] args) {

		JVector[][] axes = JVector.axes111U_ZXY;
		
		String[] lbl = new String[] {"z", "x", "y"};
		for(int i = 0; i < axes.length; i++) {
			for(int j = 0; j < axes[0].length; j++) {
				System.out.println(i + "\t" + lbl[j] + "\t" + axes[i][j]);
			}
			System.out.println("Angle between: " + axes[i][0] + " and " + axes[i][1] + "\t: " + JVector.angle(axes[i][0], axes[i][1]));
			System.out.println("Angle between: " + axes[i][0] + " and " + axes[i][2] + "\t: " + JVector.angle(axes[i][0], axes[i][2]));
			System.out.println("Angle between: " + axes[i][1] + " and " + axes[i][2] + "\t: " + JVector.angle(axes[i][1], axes[i][2]));
		}
	}
}
