package basisSets;

import io.StringConverter;
import chemistry.JAtomTools;
import defaultPackage.HookePotential;
import defaultPackage.IdealTetrahedron;
import defaultPackage.JAtom;
import defaultPackage.JVector;

public class IdealizeTetrahedron {

	private double tolerance = 0.001;
	private double bondLen = 1.91;
	private double bondAngle = 109.47;
	private double ligandDist = 2 * bondLen * Math.sin(bondAngle/2 * Math.PI / 180);
	private double ligandLigand= .1;
	private double centralLigand = 1;
	private double timeStep = 10;
	private double time = 0;
	private JAtom c;
	private JAtom[] br;
	private HookePotential brbr, cbr;
	private IdealTetrahedron target;
	boolean printForDebug = false;
	boolean printFinal = false;
	
	public IdealTetrahedron idealize() {
		time = 0;
		IdealTetrahedron idealized = (IdealTetrahedron) target.clone();
		
		br = idealized.getLigands();
		c = idealized.getCenter();

		idealizeBonds();
		idealizeAngles();
		
		return new IdealTetrahedron(c, br);
	}
	
	private void idealizeBonds() {
		for(int i = 0; i < br.length; i++) {
			JVector vec = JVector.subtract(br[i].getPosition(), c.getPosition());
			vec = JVector.multiply(vec.unit(), bondLen);
			br[i].setPosition(JVector.add(vec, c.getPosition()));
		}
	}
	
	private void idealizeAngles() {
		double[] angles = calcAngles();
		
		double stddev = stddev(angles);
		brbr = new HookePotential(35, 35, ligandDist, ligandLigand);
		cbr = new HookePotential(6, 35, bondLen, centralLigand);
		int stepNum = 0;
		int printEvery = 100;
		int numSamePotentials = 0;
		double prevPotential = calcPotential();
		double curPotential;
		double fractionDiff;
		while(numSamePotentials < 100) {
			step();
			stddev = stddev(calcAngles());
			curPotential = calcPotential();
			fractionDiff = Math.abs((prevPotential - curPotential) / curPotential);
			if(fractionDiff < tolerance) {
				numSamePotentials++;
			}
			prevPotential = curPotential;
			if(stepNum % printEvery == 0 && printForDebug) {
				System.out.println("Current potential at time step: " + time + " is: " + calcPotential());
				System.out.println("\tAngles:\t" + StringConverter.arrayToTabString(calcAngles()));
				System.out.println("\tDistances:\t" + StringConverter.arrayToTabString(calcDistances()));
			}
		}
		if(printFinal) {
			System.out.println("Current potential at time step: " + time + " is: " + calcPotential());
			System.out.println("\tAngles:\t" + StringConverter.arrayToTabString(calcAngles()));
			System.out.println("\tDistances:\t" + StringConverter.arrayToTabString(calcDistances()));
		}
	}
	
	private void step() {
		JVector[] forces = calcForces();
		move(forces);
		idealizeBonds();
		time += timeStep;
	}
	
	private void move(JVector[] forces) {
		for(int i = 0; i < forces.length; i++) {
			JVector curPosition = br[i].getPosition();
			double mass = JAtomTools.getMass(br[i].getZ());
			
			JVector move = JVector.multiply(forces[i], timeStep * timeStep / mass);
			
			br[i].setNewPos(JVector.add(curPosition, move));
		}
	}
	private JVector[] calcForces() {
		JVector[] forces = new JVector[br.length];
		for(int i = 0; i < br.length; i++) {
			forces[i] = new JVector();
			JVector cbrVec = JVector.subtract(br[i].getPosition(), c.getPosition());
			forces[i].add_inPlace(JVector.multiply(cbrVec.unit(), cbr.calcF(cbrVec.length())));
			for(int j = 0; j < br.length; j++) {
				if(i == j)
					continue;
				JVector brbrVec = JVector.subtract(br[i].getPosition(), br[j].getPosition());
				forces[i].add_inPlace(JVector.multiply(brbrVec.unit(), brbr.calcF(brbrVec.length())));
			}
		}
		
		return forces;
	}
	
	private double calcPotential() {
		double potential = 0;
		for(int i = 0; i < br.length; i++) {
			JVector cbrVec = JVector.subtract(br[i].getPosition(), c.getPosition());
			potential += cbr.calcU(cbrVec.length());
			for(int j = 0; j < br.length; j++) {
				if(i == j)
					continue;
				JVector brbrVec = JVector.subtract(br[i].getPosition(), br[j].getPosition());
				potential += brbr.calcU(brbrVec.length());
			}
		}
		return potential;
	}
	
	private double[] calcDistances() {
		double[] distances = new double[br.length];
		for(int i = 0; i < br.length; i++)
			distances[i] = JVector.subtract(br[i].getPosition(), c.getPosition()).length();
		return distances;
	}
	private double[] calcAngles() {
		int numAnglesInTetrahedron = 6;
		double[] angles = new double[numAnglesInTetrahedron];
		int idx = 0;
		JAtom br1, br2;
		for(int i = 0; i < br.length; i++) {
			br1 = br[i];
			for(int j = i+1; j < br.length; j++) {
				br2 = br[j];
				angles[idx++] = JVector.angle(JVector.subtract(br1.getPosition(), c.getPosition()), 
						JVector.subtract(br2.getPosition(), c.getPosition()));
			}
		}
		
		return angles;
	}
	/*******************/
	/** MATH FUNCTIONS */
	/*******************/
	
	private double stddev(double[] vals) {
		double avg = avg(vals);
		
		double squaredResiduals = 0;
		
		for(int i = 0; i < vals.length; i++)
			squaredResiduals = Math.pow(vals[i] - avg, 2);
		
		return Math.sqrt(squaredResiduals / vals.length);
	}
	private double avg(double[] vals) {
		double sum = sum(vals);
		return sum / vals.length;
	}
	
	private double sum(double[] vals) {
		double sum = 0;
		for(int i = 0; i < vals.length; i++)
			sum += vals[i];
		return sum;
	}
	public IdealTetrahedron getTarget() {
		return target;
	}
	public void setTarget(IdealTetrahedron target) {
		this.target = target;
	}
}
