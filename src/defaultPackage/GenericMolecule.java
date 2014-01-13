package defaultPackage;

import java.util.Vector;

public class GenericMolecule extends Molecule {

	private Vector<JAtom> atoms;
	
	public GenericMolecule(Vector<JAtom> atoms) {
		this.atoms = atoms;
	}
	
	
	@Override
	public void rotate(JVector axis, JVector origin, double phi) {
		for(JAtom atom : atoms) {
			Quaternion temp = new Quaternion(atom.getPosition());
			temp = Quaternion.rotate(temp, axis, origin, phi);
			atom.setPosition(temp.position);
		}
	}

	@Override
	public void translate(JVector amount) {
		for(JAtom atom : atoms) {
			atom.getPosition().add_inPlace(amount);
		}		
	}

	@Override
	public String toStringForAtoms(double a, double b, double c, int number,
			int aUnits, int bUnits, int cUnits) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public JAtom[] getLigands() {
		return getAtoms();
	}

	public JVector getCenterOfMass() {
		JVector weightedPosition = new JVector();
		double mass = 0;
		for(JAtom atom : atoms) {
			weightedPosition.add_inPlace(JVector.multiply(atom.getPosition(), atom.Z));
			mass += atom.Z;
		}
		weightedPosition.multiply(1./mass);
		
		return weightedPosition;
	}

	@Override
	public void setLigands(JAtom[] ligands) {
		
	}

	@Override
	public JVector[] getVectors() {
		Vector<JVector> vectors = new Vector<JVector>();
		for(JAtom atom : atoms)
			vectors.add((JVector) atom.getPosition().clone());
		
		return vectors.toArray(new JVector[vectors.size()]);
	}

	@Override
	public Object clone() {
		Vector<JAtom> cloned = new Vector<JAtom>();
		for(JAtom atom : atoms)
			cloned.add((JAtom) atom.clone());
		
		return new GenericMolecule(cloned);
	}

	@Override
	public JAtom[] getAtoms() {
		return atoms.toArray(new JAtom[atoms.size()]);
	}


	@Override
	public JAtom getCenter() {
		return atoms.get(0);
	}

	@Override
	public void moveTo(JVector newPos) {
		JAtom center = atoms.firstElement();
		JVector translationVec = JVector.subtract(newPos, center.getPosition());
		
		for(JAtom atom : atoms) {
			atom.translate(translationVec);
		}
	}
}
