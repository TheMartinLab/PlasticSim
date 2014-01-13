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
import java.io.Serializable;

/**
 * A class that defines a complex number of the form a + bi in terms of two doubles, 're' and 'im'
 * @author Eric Dill
 * 			eddill@ncsu.edu
 * @version 1.0
 *
 */
public class JComplex implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1740017228932745428L;
	
	/* Instance variables to represent the real and imaginary parts of a complex number: a + b*i = re + im * i */
	protected double re, im;

	/**
	 * Constructor, initialize the complex number to the values passed in
	 * @param re	The real part of the complex number, a
	 * @param im	The imaginary part of the complex number, b
	 */
	public JComplex(double re, double im)
	{
		this.re = re;
		this.im = im;
	}
	
	/* --------------------------------------------------------------- */
	
						/* MATHEMATICAL OPERATIONS */

	/* --------------------------------------------------------------- */
	
	/**
	 * Method to multiply two complex numbers together
	 * @param toMultiply	A complex number to be multiplied with the calling complex number
	 * @return	The multiplication of two complex numbers.  New JComplex object.
	 */
	public static JComplex multiply(JComplex one, JComplex two)
	{
		double real = one.re * two.re - one.im * two.im;
		double imaginary = one.re * two.im + one.im * two.re;
		
		return new JComplex(real, imaginary);
	}
	
	/**
	 * Method to multiply a complex number by a real number
	 * @param toMultiply	A real number
	 * @return	The calling complex number times the passed real number.  New JComplex object.
	 */
	public static JComplex multiply(JComplex one, double c) { return new JComplex(one.re * c, one.im * c); }
	
	/**
	 * Method to add two complex numbers
	 * @param toAdd	The second complex number to add
	 * @return	The sum of the calling complex number and the passed complex number.   New JComplex object.
	 */
	public static JComplex add(JComplex one, JComplex two) { return new JComplex(one.re + two.re, one.im + two.im); }
	
	/**
	 * Method to add a double to the calling complex number.
	 * @param toAdd	The double to add to the calling complex.
	 * @return	(a + ib) + (c) = new JComplex(a + c, b).  New JComplex object.
	 */
	public static JComplex add(JComplex one, double c) { return new JComplex(one.re + c, one.im); }
	
	/**
	 * Method to calculate the exponential value of this complex number in the form exp(a + bi) = exp(a)*(cos(b) + i*sin(b))
	 * @return	The complex value of exp(a + bi).   New JComplex object.
	 */
	public static JComplex exponential(JComplex one)
	{
		double ea = Math.exp(one.re);
		double real = ea * Math.cos(one.im);
		double imaginary = ea * Math.sin(one.im);
		
		return new JComplex(real, imaginary);
	}
	
	/* --------------------------------------------------------- */
	
						/* HELPER OPERATIONS */
	
	/* --------------------------------------------------------- */
	
	/**
	 * Method to calculate the conjugate of this complex number in the form conjugate(a + bi) = a - bi
	 * @return		The conjugate of this complex number.   New JComplex object.
	 */
	public JComplex conjugate() { return new JComplex(re, -1 * im); }
	
	/**
	 * Method to calculate the modulus of a complex number in the form modulus(a + bi) = (a^2 + b^2)^.5
	 * @return		The value of the modulus of this complex number
	 */
	public Double modulus() { return Math.sqrt(re * re + im * im); }
	
	
	/* ----------------------------------------------------------------- */
	
						/* GETTER AND SETTER METHODS */
	
	/* ----------------------------------------------------------------- */
	
	/**
	 * Getter method for the real part of the complex number
	 * @return	The real part of the imaginary number
	 */
	public double getRe() { return re; }
	
	/**
	 * Getter method for the imaginary part of the complex number
	 * @return	The imaginary part of the complex number
	 */
	public double getIm() { return im; }
	
	/**
	 * Setter method for the real part of the complex number
	 * @param re	The real part of the complex number
	 */
	public void setRe(double re) { this.re = re; }
	
	/**
	 * Setter method for the imaginary part of the complex number
	 * @param im	The imagainary part of the complex number
	 */
	public void setIm(double im) { this.im = im; }
	public Object clone() { return new JComplex(re, im); }
	/**
	 * toString() override to output the complex number in the form "a + bi"
	 */
	public String toString() { return re + " + " + im + "i"; }
}
