package basisSets;

import io.MyPrintStream;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Random;
import java.util.Scanner;

import defaultPackage.ComplexScatteringFactor;
import defaultPackage.DoubleLinkedListCBr4;
import defaultPackage.DoubleLinkedListVector;
import defaultPackage.IdealTetrahedron;
import defaultPackage.JAtom;
import defaultPackage.JVector;
import defaultPackage.Lattice;
import defaultPackage.LennardJonesLookup;
import defaultPackage.LennardJonesPotential;
import defaultPackage.TwoNodeCBr4;
import defaultPackage.TwoNodeVector;

public class FirstShell {

	protected int aUnits;
	protected double aSmall, aBig, aPlastic;
	protected JAtom[] carbons, bromines, box, carbons_32;
	protected IdealTetrahedron[] pseudoLattice, bigLattice,six=Lattice.getSix4barOrientations();
	protected IdealTetrahedron[][] pseudoLattices, miniLattices, fourCubes, fourBigCubes, lotsOLittleCubes, 
				rotatedBigCubes, surroundings_32, eight_options, seven_rings, sortedEightOptions, rotated_eight;
	protected LennardJonesLookup ljl;
	protected LennardJonesPotential ljp;
	protected ComplexScatteringFactor csf;
	protected IdealTetrahedron[] fcc;
	static JVector[] firstShell = JVector.firstShell;
	private File outputDirectory;
	public void readFile(File aFile) throws FileNotFoundException
	{
		FileInputStream fis = new FileInputStream(aFile);
		
		Scanner s = new Scanner(fis);

		carbons_32 = new JAtom[32];

		int cIndex = 0;
		int brIndex = 0;
		int snIndex = 0;
		int c_32 = 0;
		
		
		// figure out how many bromines and carbons I have
		while(s.hasNextLine())
		{
			String[] next = s.nextLine().split("\t");
			int Z = Integer.valueOf(next[2]);
			if(Z == 6)
				cIndex++;
			if(Z == 35)
				brIndex++;
			if(Z == 50)
				snIndex++;
		}
		carbons = new JAtom[cIndex*27];
		bromines = new JAtom[brIndex*27];
		box = new JAtom[snIndex];
		
		cIndex = 0;
		snIndex = 0;
		brIndex = 0;
		
		//double a = 21.43+ Math.sin(20.88*Math.PI/180)*21.02;
		//double b = 12.12;
		//double c = 21.02*(Math.cos(20.88*Math.PI/180));
		double a = 21.43;
		double b = 12.12;
		double c = 21.2;
		// reset the file input stream and the scanner
		fis = new FileInputStream(aFile);
		s = new Scanner(fis);
		//System.out.println("CELL\t" + a + "\t" + b + "\t" + c + "\t" + 90 + "\t" + 90 + "\t" + 90);
		//System.out.println("SPGP\tP1");
		
		double boxEdge = 17.4;
		int it_32 = 1;
		
		while(s.hasNextLine())
		{
			if(snIndex == 10)
				System.out.print("");
			a = 21.43;
			b = 12.12;
			c = 21.02;
			
			String[] next = s.nextLine().split("\t");
			int Z = Integer.valueOf(next[2]);
			
			double x = Double.valueOf(next[3]);
			double y = Double.valueOf(next[4]);
			double z = Double.valueOf(next[5]);
			
			JVector[] positions = new JVector[27];
			
			JVector pos = new JVector(x, y, z);

			int index = 0;
			
			for(int i = -1; i <= 1; i++)
			{
				for(int j = -1; j <= 1; j++)
				{
					for(int k = -1; k <= 1; k++)
					{
						pos = JVector.add(new JVector(x, y, z), new JVector(i, j, k));
						// rotate the a axis
						double xNew = JVector.dot(pos, new JVector(Math.cos(20.88*Math.PI/180), 0, 0)) ;
						double yNew = JVector.dot(pos, new JVector(0, 1, 0)) ;
						double zNew = JVector.dot(pos, new JVector(-Math.sin(20.88*Math.PI/180), 0, 1));
						// convert to cartesian coordinates
						xNew *= a;
						yNew *= b;
						zNew *= c;
						
						if(Z == 50 && (i != 0 || j != 0 || k != 0))
							continue;
						
						positions[index] = new JVector(xNew, yNew, zNew);
						
						if(i == 0 && j == 0 && k == 0 && Z == 6)
						{
							//carbons_32[c_32] = new JAtom(c_32, positions[index], csf.getAbbreviation(c_32+1));
							carbons_32[c_32] = new JAtom(6, positions[index]);
							c_32++;
						}
						index++;
					}
				}
			}

			
			//JVector newPos2 = JVector.add(newPos, new JVector(a, 0, 0));
			//JVector newPos3 = JVector.add(newPos, new JVector(a, 0, c));
			//JVector newPos4 = JVector.add(newPos, new JVector(0, 0, c));
			
			// convert to new fractional coordinates
			
			a = 1;
			b = 1;
			c = 1;
			
			if(Z == 6)
			{
				for(int i = 0; i < positions.length; i++)
				{
					//carbons[cIndex] = new JAtom(it_32, positions[i], csf.getAbbreviation(it_32));
					carbons[cIndex] = new JAtom(6, positions[i]);
					cIndex++;
				}
			}
			
			if(Z == 35)
			{
				for(int i = 0; i < positions.length; i++)
				{			
					bromines[brIndex] = new JAtom(Z, positions[i]);
					brIndex++;
				}
			}
			it_32++;
		}
		
		
		a = 21.43;
		b = 12.12;
		c = 21.02;
	}
	
	public void makeTetrahedra()
	{
		int latticePos = 0;
		pseudoLattice = new IdealTetrahedron[864];
		
		// loop through the carbons
		for(int i = 0; i < carbons.length; i++)
		{
			
			JAtom[] ligands = new JAtom[4];
			int ligandPos = 0;
			if(carbons[i] == null)
				continue;
			// loop through the bromines
			for(int j = 0; j < bromines.length; j++)
			{
				double distance = 3;
				// calculate the distance from the bromine to the carbon
				if(bromines[j] == null)
					continue;
				distance = JVector.subtract(carbons[i].getPosition(), bromines[j].getPosition()).length();
				
				// add a bromine to the carbon tetrabromide tetrahedron
				if(distance < 2.5)
				{
					ligands[ligandPos] = bromines[j];
					ligandPos++;
					
					// if there are four ligands, move on to the next carbon center
					if(ligandPos == 4) { j = bromines.length; };
				}
			}
			// make a tetrahedron from the carbon and bromine positions
			IdealTetrahedron temp = new IdealTetrahedron(carbons[i], ligands);
			
			if(ligandPos < 4)
				continue;
			// calculate the new position for the molecule
			JVector newPos = JVector.multiply(temp.getCenter().getPosition(), aSmall/17.4);
			// subtract the old position from the new position
//			temp.moveTo(newPos);
			newPos = JVector.subtract(newPos, temp.getCenter().getPosition());
			temp.translate(newPos);
			pseudoLattice[latticePos] = temp;
			latticePos++;
		}
	}
	
	public void make_32()
	{
		surroundings_32 = new IdealTetrahedron[32][13];
		//loop through the 32 carbons in the unit cell
		for(int i = 0; i < 32; i++)
		{
			JVector c_pos = carbons_32[i].getPosition();
			
			DoubleLinkedListVector c_13 = new DoubleLinkedListVector(13);
			// loop through the rest of the tetrahedron
			for(int j = 0; j < pseudoLattice.length; j++)
			{
				if(pseudoLattice[j] == null) continue;
				double length = JVector.subtract(c_pos, pseudoLattice[j].getCenter().getPosition()).length();
				c_13.insert(new TwoNodeVector(length, j, null, null));
			}
			surroundings_32[i][0] = (IdealTetrahedron) pseudoLattice[c_13.removeHead().getKey()].clone();
			JVector center = (JVector) surroundings_32[i][0].getCenter().getPosition().clone();
			surroundings_32[i][0].translate(JVector.multiply(center, -1));
			// set the surrounding tetrahedron
			for(int j = 1; j < 13; j++)
			{
				surroundings_32[i][j] = (IdealTetrahedron) pseudoLattice[c_13.removeHead().getKey()].clone();	
				surroundings_32[i][j].translate(JVector.multiply(center, -1));
			}
		}
	}
	
	public void rotate_eight()
	{
		rotated_eight = new IdealTetrahedron[eight_options.length*12][eight_options[0].length];
		JVector[] axes1 = {new JVector(1, 0, 0), new JVector(0, 1, 0), new JVector(0, 0, 1)};
		JVector[] axes2 = {new JVector(0, 0, 0), new JVector(1, 1, 1)};
		double[] phis1 = {90, 180, 270};
		double[] phis2 = {120, 240};

		JVector origin = new JVector(0, 0, 0);
		JVector axis=null, axis2=null;
		double phi=0, phi2=0;
		boolean secondaryRotation = false;
		for(int i = 0; i < eight_options.length; i++) {
			for(int j = 0; j < 12; j++) {
				switch(j){
				case 0:
					secondaryRotation = false;
					axis = axes2[0];
					phi = 0;
					break;
				case 1:
					secondaryRotation = false;
					axis = axes1[0];
					phi = phis1[0];
					break;
				case 2:
					secondaryRotation = false;
					axis = axes1[0];
					phi = phis1[1];
					break;
				case 3:
					secondaryRotation = false;
					axis = axes1[0];
					phi = phis1[2];
					break;
				case 4:
					secondaryRotation = false;
					axis = axes2[1];
					phi = phis2[0];
					break;
				case 5:
					secondaryRotation = true;
					axis = axes2[1];
					phi = phis2[0];
					axis2 = axes1[1];
					phi2 = phis1[0];
					break;
				case 6:
					secondaryRotation = true;
					axis = axes2[1];
					phi = phis2[0];
					axis2 = axes1[1];
					phi2 = phis1[1];
					break;
				case 7:
					secondaryRotation = true;
					axis = axes2[1];
					phi = phis2[0];
					axis2 = axes1[1];
					phi2 = phis1[2];
					break;
				case 8:
					secondaryRotation = false;
					axis = axes2[1];
					phi = phis2[1];
					break;
				case 9:
					secondaryRotation = true;
					axis = axes2[1];
					phi = phis2[1];
					axis2 = axes1[2];
					phi2 = phis1[0];
					break;
				case 10:
					secondaryRotation = true;
					axis = axes2[1];
					phi = phis2[1];
					axis2 = axes1[2];
					phi2 = phis1[1];
					break;
				case 11:
					secondaryRotation = true;
					axis = axes2[1];
					phi = phis2[1];
					axis2 = axes1[2];
					phi2 = phis1[2];
					break;
				}
				for(int k = 0; k < rotated_eight[0].length; k++) {
					rotated_eight[i*12+j][k] = (IdealTetrahedron) eight_options[i][k].clone();
					rotated_eight[i*12+j][k].rotate(axis, origin, phi);
					if(secondaryRotation) {
						rotated_eight[i*12+j][k].rotate(axis2, origin, phi2);
					}
				}
			}
		}
	}
	public void make_fcc_lattice()
	{
		int num_per_fcc = 13 * 14;
		int num_fcc = 1;
		int fcc_index = 0;
		Random r = new Random();
		fcc = new IdealTetrahedron[num_per_fcc * num_fcc];
		
		double aDim = 17.536/1.5;
		
		JVector translate = new JVector(0, 0, 0);
		JVector axis = new JVector(0, 1, 0);
		double phi = 20;
		JVector origin = new JVector(0, 0, 0);
		
		// make the fcc box
		for(int z = 0; z <= 2; z++)
		{
			translate.setK(z * aDim);
			for(int y = 0; y <= 2; y++)
			{
				translate.setJ(y * aDim);
				for(int x = 0; x <= 2; x++)
				{
					if((x + y + z) % 2 == 0)
					{
						translate.setI(x * aDim);
						//int selection = r.nextInt(32);
						int selection = 1;
						JVector toShift = JVector.multiply(surroundings_32[selection][0].getCenter().getPosition(), -1);
						
						// loop through the selected first shell
						for(int a = 0; a < 13; a++)
						{
							fcc[fcc_index] = (IdealTetrahedron) surroundings_32[selection][a].clone();
							fcc[fcc_index].rotate(axis, origin, phi);
							// shift based on the position of the center
							
							JVector total_translate = JVector.add(toShift, translate);
							fcc[fcc_index].translate(total_translate);
							
							fcc_index++;
						}
					}
				}
			}
		}
	}
	public String[] calc_e_surr(int index, double cutoff)
	{
		double total_e=0;
		int contacts = 0;
		
		// get the bromine ligands
		JAtom[] br = surroundings_32[index][0].getLigands();
		
		// loop through the central bromines
		for(int i = 0; i < 4; i++)
		{
			JVector br_pos = br[i].getPosition();
			// loop through the other 12 molecules
			for(int j = 1; j < 13; j++)
			{
				// get the bromines on the surroundings
				JAtom[] br_surr = surroundings_32[index][j].getLigands();
				
				// loop through the ligand bromines
				for(int k = 0; k < 4; k++)
				{
					double len = JVector.subtract(br_pos, br_surr[k].getPosition()).length();
					if(len < cutoff)
					{
						total_e += ljp.calcU(len);
						//total_e += ljl.lookupPotential(len);
						contacts++;
					}
						
				}
			}
		}
		String[] temp = {"" + total_e/2, "" + contacts};
		
		return temp;
	}
	public String makeAtomsHeader(double a)
	{
		String header = "";
		
		//double a = aUnits * aSmall/2;
		
		header += "CELL " + a + " " + a + " " + a + " " + 90 + " " + 90 + " " + 90 + "\n";
		header += "SPGP P1\n";
		
		return header;
	}
	
	public void print_fcc()
	{
		String fname = "fcc.inp";
		FileOutputStream fos = null;
		PrintStream ps = null;

		String output_name = fname;
		try {
			fos = new FileOutputStream(output_name);
			ps = new PrintStream(fos);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		ps.println(makeAtomsHeader(1));
		
		for(int i = 0; i < fcc.length; i++)
		{
			ps.println(fcc[i].toStringForAtoms(1, 1, 1, 1, 1, 1, 1));
		}
	}
	public void print_32() throws IOException
	{
		String fname = "_of_32_tetrahedron.inp";
		String name = "";
		FileOutputStream fos = null;
		PrintStream ps = null;

		for(int i = 0; i < 32; i+=4)
		{
			name = "" + i + fname;
			fos = new FileOutputStream(name);
			ps = new PrintStream(fos);
			ps.println(makeAtomsHeader(1));
			for(int j = 0; j < 13; j++)
			{/*
				ps.println("6," + surroundings_32[i][j].getCenter());
				for(int k = 0; k < 4; k++)
				{
					JAtom[] ligands = surroundings_32[i][j].getLigands();
					ps.println("35," + ligands[k].getPosition());
				}
				*/
				ps.println(surroundings_32[i][j].toStringForAtoms(1, 1, 1, 1, 1, 1, 1));
			}
		}
	}
	
	public void align()
	{
		eight_options = new IdealTetrahedron[8][13];
		JVector origin = new JVector(0, 0, 0);
		JVector axis = new JVector(0, 1, 0);
		double phi = 35;
		int index = 0;
		for(int i = 0; i < 32; i+=4)
		{
			for(int j = 0; j < 13; j++)
			{
				surroundings_32[i][j].rotate(axis, origin, phi);
			}
		}
		axis = new JVector(0, 0, 1);
		phi = 45;
		for(int i = 0; i < 32; i+=4)
		{
			for(int j = 0; j < 13; j++)
			{
				surroundings_32[i][j].rotate(axis, origin, phi);
				eight_options[index][j] = (IdealTetrahedron) surroundings_32[i][j].clone();
			}
			index++;
		}
	}
	public void sortEight(IdealTetrahedron[][] lattice) {
		DoubleLinkedListCBr4 dll = new DoubleLinkedListCBr4(5);
		sortedEightOptions = new IdealTetrahedron[lattice.length][lattice[0].length];
		double dist;
		for(int i = 0; i < lattice.length; i++) {
			for(int k = 0; k < firstShell.length; k++) {
				dll.clear();
				for(int j = 0; j < lattice[0].length; j++) {
					dist = JVector.subtract(lattice[i][j].getCenter().getPosition(), JVector.multiply(firstShell[k], 8.82/2)).length();
					dist = Math.abs(dist);
					dll.insert(new TwoNodeCBr4(lattice[i][j], dist, null, null));
				}
				sortedEightOptions[i][k] = (IdealTetrahedron) dll.removeHead().getValue().clone();
			}
		}
	}
	public void get_7_rings()
	{
		seven_rings = new IdealTetrahedron[8][7];
		
		// get the first two from each second shell
		for(int i = 0; i < 8; i++)
			for(int j = 0; j < 2; j++)
				seven_rings[i][j] = (IdealTetrahedron) eight_options[i][j].clone();

		DoubleLinkedListVector list = new DoubleLinkedListVector(11);
		JVector cross, pos_1, pos_2, pos_test;
		double val;
		int which_mol;
		
		// get the other 5
		for(int i = 0; i < 8; i++)
		{
			list = new DoubleLinkedListVector(11);
			pos_1 = JVector.subtract(seven_rings[i][1].getCenter().getPosition(), seven_rings[i][0].getCenter().getPosition());
			
			// find one on the same ring
			for(int j = 2; j < 13; j++)
			{
				pos_test = JVector.subtract(eight_options[i][j].getCenter().getPosition(), seven_rings[i][0].getCenter().getPosition());
				val = Math.abs(Math.acos(JVector.dot(pos_1, pos_test) / (pos_1.length() * pos_test.length())));
				val *= 180/Math.PI;
				list.insert(new TwoNodeVector(val, j, null, null));
			}
			which_mol = list.removeHead().getKey();
			System.out.println("third: " + which_mol);
			seven_rings[i][2] = (IdealTetrahedron) eight_options[i][which_mol].clone();
			pos_2 = JVector.subtract(eight_options[i][which_mol].getCenter().getPosition(), eight_options[i][0].getCenter().getPosition());
			
			cross = JVector.cross(pos_1, pos_2);
			list = new DoubleLinkedListVector(10);
			
			for(int j = 2; j < 13; j++)
			{
				if(j == which_mol)
					continue;
				pos_test = JVector.subtract(eight_options[i][j].getCenter().getPosition(), seven_rings[i][0].getCenter().getPosition());
				val = Math.abs(JVector.dot(cross, pos_test));
				System.out.println("j: " + j + "\tval: " + val);
				list.insert(new TwoNodeVector(val, j, null, null));
			}
			// fill in the remaining 4
			for(int j = 3; j < 7; j++)
			{
				which_mol = list.removeHead().getKey();
				System.out.println(which_mol);
				seven_rings[i][j] = (IdealTetrahedron) (eight_options[i][which_mol]).clone();
			}
		}
	}
	
	public void displace(double amount) {
		int whichOrientation;
		DoubleLinkedListVector dllv1 = new DoubleLinkedListVector(10);
		DoubleLinkedListVector dllv2 = new DoubleLinkedListVector(10);
		IdealTetrahedron one, two;
		double angle, totalAngle;
		JVector trans = null, origin = new JVector(0, 0, 0);
		for(int i = 0; i < eight_options.length; i++) {
			dllv2.clear();
			one = eight_options[i][0];
			for(int k = 0; k < six.length; k++) {
				two = six[k];
				totalAngle = 0;
				for(int a = 0; a < one.getLigands().length; a++) {
					dllv1.clear();
					for(int b = 0; b < two.getLigands().length; b++) {
						angle = JVector.angle(one.getLigands()[a].getPosition(), two.getLigands()[b].getPosition());
						dllv1.insert(new TwoNodeVector(angle, b, null, null));
					}
					totalAngle += dllv1.removeHead().getValue();
				}
				dllv2.insert(new TwoNodeVector(totalAngle, k, null, null));
			}
			whichOrientation = dllv2.removeHead().getKey();
			switch(whichOrientation) {
			case 0:
			case 1:
				trans = (JVector) JVector.xy.clone(); 
				break;
				
			case 2:
			case 3:
				trans = (JVector) JVector.yz.clone();
				break;
				
			case 4:
			case 5:
				trans = (JVector) JVector.xz.clone();
				break;
			}
			try {
				trans = JVector.multiply(trans.unit(), amount);
			} catch (Exception e) {
			e.printStackTrace();
			}
			for(int j = 0; j < eight_options.length; j++) {
				eight_options[i][j].translate(trans);
			}
		}
		
	}
	public void mono_to_cartesian()
	{
		for(int i = 0; i < surroundings_32.length; i++) {
			File aFile = new File(outputDirectory + File.separator + "Molecule " + (i+1) + " of 32 and its surroundings.inp");
			aFile.getParentFile().mkdirs();
			MyPrintStream mps = new MyPrintStream(aFile);
			mps.println(makeAtomsHeader(1));
			
			for(int j = 0; j < surroundings_32[i].length; j++)
				mps.println(surroundings_32[i][j].toStringForAtoms(1, 1, 1, 1, 1, 1, 1));
			
			mps.close();
		}
	}
	public void print(IdealTetrahedron[][] lattice, String aFile)
	{
		FileOutputStream fos = null;
		PrintStream ps = null;
		try {
			fos = new FileOutputStream(aFile);
			ps = new PrintStream(fos);
		} catch (IOException e) {
			e.printStackTrace();
		}
		ps.println(lattice.length*lattice[0].length);
		ps.println(lattice.length);
		ps.println(lattice[0].length);
		for(int i = 0; i < lattice.length; i++) {
			ps.println(i);
			for(int j = 0; j < lattice[0].length; j++) {
				ps.print(lattice[i][j].toStringForXYZ());
			}
			
		}
	}
	public static void main(String[] args)
	{
		FirstShell m = new FirstShell();
		
		// deltaH_vap = 
		m.ljp = new LennardJonesPotential(35, 35, 3.3854, 0.029809248);
		
		m.ljl = new LennardJonesLookup(.001, m.ljp);
		
		m.csf = new ComplexScatteringFactor();
		
		m.aUnits = 2;
		m.aSmall = 2*8.82;
		m.aBig = 22;
		m.aPlastic = 8.82;
		
		
		m.outputDirectory = new File("output");
		
		m.outputDirectory = m.outputDirectory.getAbsoluteFile();
		
		try {

			m.readFile(new File("monoclinic crystallographic coords.txt")); 
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}

		m.makeTetrahedra();
		
		m.make_32();
		IdealTetrahedron[][] m32 = m.surroundings_32;
		m.mono_to_cartesian();
		FileOutputStream fos = null;
		PrintStream ps = null;
		String output_name = "ener_surroundings.txt";
		try {
			fos = new FileOutputStream(output_name);
			ps = new PrintStream(fos);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		// print the contacts first
		String[] temp = null;
		for(int j = 0; j < 20; j++)
		{
			ps.print(j);
			for(int i = 0; i < 32; i++)
			{
				temp = m.calc_e_surr(i, j);
				ps.print("\t" + temp[1]);
			}
				
			ps.println();
		}
		
		ps.println("energy:");
		
		// print the total energy
		for(int j = 0; j < 20; j++)
		{
			ps.print(j);
			for(int i = 0; i < 32; i++)
			{
				temp = m.calc_e_surr(i, j);
				ps.print("\t" + temp[0]);
			}
				
			ps.println();
		}
		m.align();
		m.displace(.4);
		m.rotate_eight();
		// sort the first shell molecules such that they are in the order specified by the JVector.firstShell Vector order
		m.sortEight(m.rotated_eight);
		
		try {
			m.print_32();
			m.print(m.sortedEightOptions, "sorted_eight.lattice");
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		m.make_fcc_lattice();
		m.print_fcc();
		m.print(m.eight_options, "eight_options.lattice");
	}
}
