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
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashMap;

public class LennardJonesLookup {

	private HashMap<Double, Double> energyLookup;
	
	private HashMap<Double, Double> forceLookup;
	
	private LennardJonesPotential ljPotential;
	
	private double precision;
	
	public LennardJonesPotential getLennardJones() { return ljPotential; }
	
	public LennardJonesLookup(double precision, LennardJonesPotential ljPotential)
	{
		this.ljPotential = ljPotential;
		
		this.precision = precision;
		
		//readPotentials();
		createPotentials();
	}

	public double lookupForce(double distance)
	{
		if(distance + precision > ljPotential.sigma * 2.5 || distance == 0)
			return 0;
		
		distance = Math.rint(distance / precision) * precision;
		
		double potential=0;
		try {
			potential = forceLookup.get(distance);
		} catch(NullPointerException npe) {
			potential = ljPotential.calcF(distance);
			forceLookup.put(distance, potential);
		}
	
		return potential;
	}
	
	public double lookupPotential(double distance)
	{
		if(distance + precision > ljPotential.sigma * 2.5 || distance == 0)
			return 0;
		
		distance = Math.rint(distance / precision) * precision;
		double potential=0;
		try {
			potential = energyLookup.get(distance);
		} catch(NullPointerException npe) {
			potential = ljPotential.calcU(distance);
			energyLookup.put(distance, potential);
		}

		return potential;
	}
	
	@SuppressWarnings("unchecked")
	public void readPotentials()
	{
		/* Declare the input streams */
		FileInputStream fis = null;
		ObjectInputStream ois = null;
		
		/* Declare the file name*/
		String fileNameE = precision + "_" + ljPotential.toString() + "_LennardJonesPotential&ForceLookup.eric";
		
		/* Find out where the path is */
		File temp = new File(".");
		
		String path = temp.getAbsolutePath();
		
		path = path.substring(0, path.length()-1) + "PotentialLookups/" + fileNameE;
		
		boolean tryAgain = false;
		
		do
		{
			/* Initialize the input streams */
			try {
				fis = new FileInputStream("Br-Br_default.eric");
				ois = new ObjectInputStream(fis);
				tryAgain = false;
			} catch (FileNotFoundException e) {
				writePotentials();
				tryAgain = true;
			} catch (IOException e) {
				System.out.println("Could not create an ObjectInputStream");
				e.printStackTrace();
			} 
		} while(tryAgain == true);

		/* input the objects */
		do
		{
			try {
				energyLookup = (HashMap<Double, Double>) ois.readObject();
				forceLookup = (HashMap<Double, Double>) ois.readObject();
				tryAgain = false;
			} catch (IOException e) {
				writePotentials();
				readPotentials();
				tryAgain = true;
			} catch (ClassNotFoundException e) {
				e.printStackTrace();
			}
		} while(tryAgain == true);

	}
	public void createPotentials() {
		int numDataPoints = (int) (2.5 * ljPotential.sigma / precision);
		
		double maxValue = numDataPoints * precision;
		
		energyLookup = new HashMap<Double, Double>(numDataPoints);
		
		forceLookup = new HashMap<Double, Double>(numDataPoints);
		
		for(double i = 0; i <= maxValue; i += precision)
		{
			double distance = Math.rint(i / precision) * precision;
			
			double potential = ljPotential.calcU(distance);
			
			double force = ljPotential.calcF(distance);
			
			energyLookup.put(distance, potential);
			
			forceLookup.put(distance, force);
		}
	}
	/**
	 * 
	 * @param precision
	 */
	public void writePotentials()
	{
		int numDataPoints = (int) (2.5 * ljPotential.sigma / precision);
		
		double maxValue = numDataPoints * precision;
		
		energyLookup = new HashMap<Double, Double>(numDataPoints);
		
		forceLookup = new HashMap<Double, Double>(numDataPoints);
		
		for(double i = 0; i <= maxValue; i += precision)
		{
			double distance = Math.rint(i / precision) * precision;
			
			double potential = ljPotential.calcU(distance);
			
			double force = ljPotential.calcF(distance);
			
			energyLookup.put(distance, potential);
			
			forceLookup.put(distance, force);
		}
		
		/* Declare the output stream */
		FileOutputStream fos = null;
		ObjectOutputStream oos = null;
		
		/* Declare the file name*/
		String fileName = precision + "_" + ljPotential.toString() + "_LennardJonesPotential&ForceLookup.eric";
		
		/* Find out where the path is */
		File temp = new File(".");
		
		String path = temp.getAbsolutePath();
		
		/* Create the path to the file */
		path = path.substring(0, path.length()-1) + "/PotentialLookups/" + fileName;

		/* initialize the output stream */
		try {
			fos = new FileOutputStream("Br-Br_default.eric");
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
			oos.writeObject(energyLookup);
			oos.writeObject(forceLookup);
		} catch (IOException e) {
			System.out.println("Cannot write the objects to disk");
			e.printStackTrace();
		}

		
	}
	
	public double getPrecision() { return precision; }
}
