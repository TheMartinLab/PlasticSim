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
 * This class defines a quaternion in the form of (s, i, j, k).  Available operations: Rotation, Normalize, 
 * Multiply by a constant, Multiply two Quats, Determine Quaternion conjugate, Add 2 quaternions, Subtract 2 quaternions
 * 
 * @author Eric D. Dill eddill@ncsu.edu
 * @author James D. Martin jdmartin@ncsu.edu
 * Copyright © 2010-2013 North Carolina State University. All rights reserved
 */

public class Quaternion implements Cloneable
{

	/* INSTANCE VARIABLES */
	protected double s;

	protected JVector position;
	
	public Quaternion(double s, JVector v)
	{
		position = (JVector) v.clone();
		
		this.s = s;
	}
	
	public Quaternion(JVector v)
	{
		position = (JVector) v.clone();
		
		s = 0;
	}
	
	public Quaternion()
	{
		position = new JVector();
		
		s = 0;
	}
	

	/* --------------------------------------------------------------- */

						/* CLASS OPERATIONS */

	/* --------------------------------------------------------------- */

	/**
	 * This method calculates the clockwise rotation of a quaternion respective to a JVector origin
	 * relative to a JVector axis by an angle phi.
	 * @param axis		The JVector representing the axis around which the calling quaternion is to be rotated.
	 * @param origin	The JVector representing the origin.
	 * @param phi		The angle which the calling quaternion is to be rotated.
	 * @return 			The rotated quaternion. New Quaternion object.
	 */
	public static Quaternion rotate(Quaternion toRotate, JVector vectorAxis, JVector vectorOrigin, double phi) 
	{
		/* Create a quaternion representation of the axis */
		Quaternion axis = new Quaternion(vectorAxis);
		
		/* Create a quaternion representation of the origin */
		Quaternion origin = new Quaternion(vectorOrigin);
		
		/* Create a temporary quaternion variable */
		Quaternion temp = new Quaternion();
		
		/* Calculate the cosine of the angle phi */
		double cosPhi = Math.cos(phi * Math.PI / 360);
		
		/* Calculate the sine of the angle phi */
		double sinPhi = Math.sin(phi * Math.PI / 360);
		
		/* normalize the axis */
		try {
			axis = axis.unit();
		} catch (Exception e) {
			return toRotate;
		}
		
		/* turn the axis into the rotation quaternion */
		axis = multiply(axis, sinPhi);
		axis.s = cosPhi;
		
		Quaternion p = subtract(toRotate, origin);
		
		temp = cross(p, axis.conjugate());
		
		temp = cross(axis, temp);
		
		/* Put the quaternion back where it started */
		temp = add(temp, origin);
		
		/* Return the newly rotated quaternion */
		return temp;
	}
	
	public static Quaternion rotate(Quaternion toRotate, JVector vectorAxis, JVector vectorOrigin, double phi, Quaternion[] space) 
	{
		/* Create a quaternion representation of the axis */
		space[0] = new Quaternion(vectorAxis);
		
		/* Create a quaternion representation of the origin */
		space[1] = new Quaternion(vectorOrigin);
		
		/* Create a temporary quaternion variable */
		space[2] = new Quaternion();
		
		/* Calculate the cosine of the angle phi */
		double cosPhi = Math.cos(phi * Math.PI / 360);
		
		/* Calculate the sine of the angle phi */
		double sinPhi = Math.sin(phi * Math.PI / 360);
		
		/* normalize the axis */
		try {
			space[0] = space[0].unit();
		} catch (Exception e) {
			return toRotate;
		}
		
		/* turn the axis into the rotation quaternion */
		space[0] = multiply(space[0], sinPhi);
		space[0].s = cosPhi;
		
		space[3] = subtract(toRotate, space[1]);
		
		space[2] = cross(space[3], space[0].conjugate());
		
		space[2] = cross(space[0], space[2]);
		
		/* Put the quaternion back where it started */
		space[2] = add(space[2], space[1]);
		
		/* Return the newly rotated quaternion */
		return space[2];
	}
	
	
	/**
	 * Method to multiply a Quaternion by a real
	 * @param q A quaternion
	 * @param c A constant
	 * @return	A new quaternion object equal to q * c
	 */
	public static Quaternion multiply(Quaternion q, double c)
	{
		return new Quaternion(q.s * c, JVector.multiply(q.position, c));
	}
	
	/**
	 * Method to multiply two Quaternions together and return the product
	 * @param q1 A quaternion
	 * @param q2 Another quaternion
	 * @return	A new quaternion object equal to q1 x q2
	 */
	public static Quaternion cross(Quaternion q1, Quaternion q2)
	{
		double a = q1.s * q2.s - q1.position.i * q2.position.i - q1.position.j * q2.position.j - q1.position.k * q2.position.k;
		
		double b = q1.s * q2.position.i + q1.position.i * q2.s + q1.position.j * q2.position.k - q1.position.k * q2.position.j;
		
		double c = q1.s * q2.position.j - q1.position.i * q2.position.k + q1.position.j * q2.s + q1.position.k * q2.position.i;
		
		double d = q1.s * q2.position.k + q1.position.i * q2.position.j - q1.position.j * q2.position.i + q1.position.k * q2.s;
		
		JVector newPosition = new JVector(b, c, d);
		
		return new Quaternion(a, newPosition);
	}
	
	/**
	 * Method to subtract two quaternions
	 * @param q1 A quaternion
	 * @param q2 Another quaternion
	 * @return	A new quaternion object equal to q1 - q2
	 */
	public static Quaternion subtract(Quaternion q1, Quaternion q2)
	{
		double newS = q1.s - q2.s;
		
		JVector newPosition = JVector.subtract(q1.position, q2.position);
		
		return new Quaternion(newS, newPosition);		
	}
	
	/**
	 * Method to add two Quaternions
	 * @param q1 A quaternion
	 * @param q2 Another quaternion
	 * @return A new quaternion object equal to q1 + q2 
	 */
	public static Quaternion add(Quaternion q1, Quaternion q2)
	{
		return new Quaternion(q1.s + q2.s, JVector.add(q1.position, q2.position));
	}
	
	
	/* --------------------------------------------------------- */
	
						/* INSTANCE OPERATIONS */
	
	/* --------------------------------------------------------- */
	
	
	/**
	 * Method to determine the conjugate of the calling quaternion.
	 * @return	A new quaternion object equal to the conjugate of the calling Quaternion object
	 */
	public Quaternion conjugate()
	{
		return new Quaternion(s, JVector.multiply(position, -1));
	}
	
	/**
	 * Method to determine the unit Quaternion of the calling Quaternion
	 * @return	The unit Quaternion of the calling Quaternion.  New Quaternion object
	 * @throws Exception An exception is thrown if the length variable is zero
	 */
	public Quaternion unit() throws Exception
	{
		double i = position.i;
		
		double j = position.j;
		
		double k = position.k;
		
		double length = Math.sqrt(s*s + i*i + j*j + k*k);
		
		if(length == 0)
			throw new Exception("The length of the Quaternion is zero");
		
		return new Quaternion(s / length, JVector.multiply(position, 1/length));
	}
	
	/* --------------------------------------------------------------- */
	
						/* GETTER & SETTER METHODS */
	
	/* --------------------------------------------------------------- */

	/**
	 * Method to update the i, j, k parameters of the Quaternion from a vector
	 * @param i	The i coordinate
	 * @param j	The j coordinate
	 * @param k	The k coordinate
	 */
	public void setLocation(JVector v)
	{
		position = v;
	}
	
	/**
	 * Get a new JVector object corresponding to the (i, j, k) position of the quaternion
	 * @return A new JVector object corresponding to the (i, j, k) position of the quaternion
	 */
	public JVector getIJK() { return (JVector) position.clone(); }
	
	/**
	 * Setter method for the s-coordinate
	 * @param s	The s-coordinate
	 */
	public void setS(double s) { this.s = s; }
	
	/**
	 * Getter method for the s-coordinate
	 * @return	The s-coordinate
	 */
	public double getS() { return s; }
	
	/**
	 * toString() override method to output the Quaternion in the form (s, i, j, k)
	 */
	@Override
	public String toString() { return + s + " (" + position.i + ", " + position.j + ", " + position.k + ")"; }
	
	public Object clone()
	{
		JVector newPosition = (JVector) position.clone();
		
		return new Quaternion(s, newPosition);
	}
}
