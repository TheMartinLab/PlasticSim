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
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashMap;
import java.util.Scanner;

public class ComplexScatteringFactor {
	
	private double[][] elementConstants;
	
	private String[][] elementNameAbbrev;
	
	public ComplexScatteringFactor() { readElementObjects(); }

	public HashMap<Double, JComplex> buildHashMap(double qMax, double qStep, double wavelength, int Z, JVector vx, JVector vy)
	{
		JComplex complexF = null;
		
		try {
			complexF = getComplexF(1, wavelength, Z);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		double mapSize = qMax / qStep / .75;
		
		HashMap<Double, JComplex> hm = new HashMap<Double, JComplex>((int) mapSize);
		
		for(double qy = 0; qy < qMax; qy+= qStep)
		{
			for(double qx = 0; qx < qMax; qx += qStep)
			{
				JVector Qx = JVector.multiply(vx, qx);
				
				JVector Qy = JVector.multiply(vy, qy);
				
				JVector Q = JVector.add(Qx, Qy);
				
				double key = Q.length();
				
				if(!hm.containsKey(key))
				{
					JComplex temp = generateF0(Q, Z, complexF);

					hm.put(key, temp);
				}
			}
		}
		return hm;
	}

	/**
	 * This method generates the complex scattering factors based on the Cromer-Mann coefficents: <br>
	 * <br> reference: Cromer, D. T.; Mann, J. B., X-RAY SCATTERING FACTORS COMPUTED FROM NUMERICAL HARTREE-FOCK WAVE FUNCTIONS. Acta Crystall a-Crys 1968, A 24, 321-&. <br>
	 * <br> website: http://www.ruppweb.org/xray/comp/scatfac.htm <br>
	 * 
	 * @param Q The reciprocal lattice vector
	 * @param Z The atomic number
	 * @param energyFormFactor The anomalous dispersion from the "getComplexF(int, double, int)" method
	 * @return The total scattering factor, f = f0 + f' + i*f"
	 */
	public JComplex generateF0(JVector Q, int Z, JComplex energyFormFactor)
	{
		JComplex f0 = new JComplex(0, 0);
		
		// http://www.ruppweb.org/xray/comp/scatfac.htm
		double fTemp = elementConstants[Z-1][9];
		double ai, bi;
		for(int i = 1; i < 5; i++)
		{
			ai = elementConstants[Z-1][i];
			bi = elementConstants[Z-1][i+4];
			
			fTemp += ai * Math.exp(-1 * bi * Math.pow((Q.length()) / (4 * Math.PI), 2));
		}
		
		f0 = JComplex.add(energyFormFactor, fTemp);
		
		return f0;
	}
	
	
	/**
	 * This method searches through the "sf" folder which contains the complex scattering factors at a number of energies in eV
	 * for the desired energy.  These files were obtained from the Berkeley National Lab.  url http://henke.lbl.gov/optical_constants/asf.html
	 * The files that will be searched only have energy values between 10eV and 30keV.  If a selection of > 30keV is made then a complex value of Z + 0*i will be returned.
	 * @param energyOrWavelength	set to 0 for an Energy in eV input.  set to 1 for a wavelength in Angstroms input
	 * @param findMe	The energy or wavelength that is to be found in the list
	 * @param Z	The number of electrons in the element
	 * @return	The complex scattering factor for the desired energy or wavelength
	 * @throws Exception	If the energyOrWavelength toggle is < 0 or > 1 then an exception is thrown.  
	 * 						If Z < 0 or > 92 an exception is thrown.  
	 * 						If E < 10 eV or lambda > 1238.174 Angstroms an exception is thrown
	 */
	public JComplex getComplexF(int energyOrWavelength, double findMe, int Z) throws Exception
	{
		/* make sure the input is 0 or 1 */
		if(energyOrWavelength < 0 || energyOrWavelength > 1)
			throw new Exception("Your selection of " + energyOrWavelength + " does not fall within the constraints." + 
					"\nIf you want energy, set the energyOrWavelength value to 0.\nIf you want wavelength, set the energyOrWavelength value to 1.");
		
		/* make sure the 0 < Z < 92 */
		if(Z < 0 || Z > 92)
			throw new Exception("Your selection of " + Z + " does not fall within the current number of elements.  Please choose a Z such that 0 < Z < 92");
		
		/* Declare the file reader */
		FileReader fr = null;
		
		/* define Planck's constant in terms of eV*s */
		double h = 4.13566733e-15;
		
		/* define the speed of light in terms of Angstroms / second */
		double c = 2.99792458e18;
		
		/* If the user passes in a wavelength (energyOrWavelength == 1), convert the wavelength to an energy in eV */
		double energyDesired;
		if(energyOrWavelength != 0)
			energyDesired = h * c / findMe; 
		else
			energyDesired = findMe;
		
		/* Check to make sure that the energy or wavelength passed in is within the parameters available to this method. */
		if(energyDesired < 10)
			throw new Exception("Your choice of an energy or wavelength of " + energyDesired + " is too low. \nPlease choose an energy such that " +
					"E > 10eV or a wavelength such that lambda < 1238.174 Angstroms."); 
		
		/* If the energy that is passed in is greater than 30keV or less than .413 Angstroms, return a real scattering value of the Z value and 0 for the complex scattering.
		 * This is a reasonable approximation.  If more exact numbers are needed, please reference http://www.nist.gov/physlab/data/ffast/index.cfm */
		if(energyDesired > 30000)
			return new JComplex(Z, 0);
		
		
		/* get the abbreviation for the element based on the index i and make it lower case */
		String abbreviation = elementNameAbbrev[Z-1][0].toLowerCase();
		
		/* find the location of the program by creating a new file "." and getting the path.  Then take the substring that includes the path less the file name */
		String findLocation = (new File(".").getAbsolutePath());
		
		findLocation = findLocation.substring(0, findLocation.length()-1);
		
		/* Initialize the file reader */
		try {
			fr = new FileReader(new File(findLocation + "/sf/" + abbreviation + ".nff"));
		} catch (FileNotFoundException e) {
			System.out.println("Cannot find the .nff files");
			e.printStackTrace();
		}
		
		/* Declare a new scanner to parse the file 'name.nff' */
		Scanner s = new Scanner(fr);
		
		/* skip the first line of the file */
		if(s.next().compareTo("E(eV)") == 0)
			s.nextLine();
		
		/* declare two variables to hold the current and previous values of the energy to compare them */ 
		double energyCurrent = 0;
		double energyPrevious = 0;
		
		/* declare two JComplex numbers to hold the current and previous values of the complex scattering factors */
		JComplex current = new JComplex(0., 0.);
		JComplex previous = new JComplex(0., 0.);
		
		/* if findMe is between two energies then take a percentage of the surrounding energies or wavelengths */
		double howMuchOf1 = 0;
		double howMuchOf2 = 0;
		
		/* declare a boolean toggle to break out of the while loop */
		boolean found = false;
		
		/* Loop through the file as long as the desired energy has not been found*/
		while(found == false)
		{
			/* set the previous energy value to what used to be the current energy value */
			energyPrevious = energyCurrent;
			
			/* set the previous scattering factor to what used to be the current scattering factor */
			previous.setRe(current.getRe());
			previous.setIm(current.getIm());
			
			/* read energy */
			energyCurrent = Double.valueOf(s.next());
			
			/* read real and imaginary parts of the scattering factor */
			current.setRe(Double.valueOf(s.next()));
			current.setIm(Double.valueOf(s.next()));
			
			/* test the current energy value to the target energy value */
			if(energyCurrent > energyDesired)
			{
				/* get out of the loop */
				found = true;
				
				/* how far away is the first energy */
				howMuchOf1 = Math.abs(energyPrevious - energyDesired) / (Math.abs(energyPrevious - energyDesired) + Math.abs(energyCurrent - energyDesired));
				
				/* how far away is the second energy */
				howMuchOf2 = Math.abs(energyCurrent - energyDesired)/ (Math.abs(energyPrevious - energyDesired) + Math.abs(energyCurrent - energyDesired));
	
			}
	
			if(found == false)
				s.nextLine();
		}
	
		
		previous = JComplex.multiply(previous, howMuchOf1);
		
		/* set the second complex scattering factor to be how far away from the second energy your selection is times 
		   the complex scattering factor at the second energy */
		current = JComplex.multiply(current, howMuchOf2);
	
		/* Test to see if the first complex scattering factor is near an absorption edge */
		if(previous.getRe() < -1 * Z)
			previous.setRe(0);
		
		/* Test to see if the second complex scattering factor is near an absorption edge */
		if(current.getRe() < -1 * Z)
			current.setRe(0);	
		
		/* return the complex scattering factor to be how far away from the bracketing energies your selection is times 
		   the complex scattering factor */
		return JComplex.add(previous, current);
	}

	/**
	 * Method to read the element abbreviation and name object and the element Z, weight, and CM coefficents from disk
	 */
	public void readElementObjects()
	{

		/* Declare the input streams */
		FileInputStream fis = null;
		ObjectInputStream ois = null;
		
		boolean tryAgain = false;
		
		do
		{
			/* Initialize the input streams */
			try {
				fis = new FileInputStream("JavaElementArrayObject");
				ois = new ObjectInputStream(fis);
				tryAgain = false;
			} catch (FileNotFoundException e) {
				writeElementInfoAsObjects();
				tryAgain = true;
			} catch (IOException e) {
				System.out.println("Could not create an ObjectInputStream");
				e.printStackTrace();
			} 
		} while(tryAgain == true);

		/* input the objects */
		try {
			elementNameAbbrev = (String[][]) ois.readObject();
			elementConstants = (double[][]) ois.readObject();
		} catch (IOException e) {
			System.out.println("Could not read the file");
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			System.out.println("Could not find the double[][] or String[][] class");
			e.printStackTrace();
		}
	}
	
	/**
	 * Generate the Element name and abbreviation object and the Element atomic weight, Z, and Cromer-Mann coefficients
	 */
	public void writeElementInfoAsObjects()
	{
		
		/* Create the String array to hold the elements */
		String[][] elements = new String[118][2];
		
		/* Create the double array to hold the atomic weight, Z, and CM Coefficients */
		double[][] elementInfo = new double[118][10];
		

		/* Initialize the file reader */
		FileReader fr = null;
		try {
			fr = new FileReader("ElementListAndFormFactors.txt");
		} catch (FileNotFoundException e) {
			System.out.println("Cannot find the ElementListAndFormFactors.txt file");
			e.printStackTrace();
		}
		
		Scanner s = new Scanner(fr);
			
		/* skip the first line of the file */
		s.nextLine();
		
		/* create an integer variable to denote the line of the file */
		int index = 0;
		
		/* Loop through the file as long as there is a next line.
		 * There are Element abbreviations and names up to Z = 104 and CM coefficents up to Z = 46
		 */
		while(s.hasNext())
		{
			/* read abbreviation */
			elements[index][0] = s.next();
			
			/* read name */
			elements[index][1] = s.next();
			
			/* read atomic weight */
			elementInfo[index][0] = Double.valueOf(s.next());
			
			/* skip the Z value since Z = index + 1 */
			s.next();
			
			/* only get the CM coefficients if Z < 47 */
			if(index < 46)
				/* loop through and input the CM coefficients */
				for(int j = 1; j < 10; j++)
					elementInfo[index][j] = Double.valueOf(s.next());
			
			s.nextLine();
			index++;
		}
		
		/* Declare the output stream */
		FileOutputStream fos = null;
		ObjectOutputStream oos = null;
		
		/* initialize the output stream */
		try {
			fos = new FileOutputStream("JavaElementArrayObject");
			oos = new ObjectOutputStream(fos);
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create the JavaElementArrayObject File");
			e.printStackTrace();
		} catch (IOException e) {
			System.out.println("Could not create an ObjectOutputStream");
			e.printStackTrace();
		}
		
		/* write the output stream */
		try {
			oos.writeObject(elements);
			oos.writeObject(elementInfo);
		} catch (IOException e) {
			System.out.println("Cannot write the objects to disk");
			e.printStackTrace();
		}

		
	}

	public String getName(int z) { return elementNameAbbrev[z-1][1]; }

	public String getAbbreviation(int z) { return elementNameAbbrev[z-1][0]; }
	
	public int getZ(String abbrev)
	{
		boolean found = false;
		
		int Z = 0;
		
		while(!found)
			if(elementNameAbbrev[Z][0].toLowerCase().compareTo(abbrev.toLowerCase()) == 0)
			{
				Z++;
				return Z;
			}
		return Integer.MAX_VALUE;
	}
}
