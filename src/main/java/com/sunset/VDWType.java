/**
 * 
 */
package com.sunset;



/**
 * @author Oleg Ursu
 *
 */
public class VDWType extends Pattern {
	
	public Bool polar;
	public double radius;
	
	public VDWType(String smarts, Bool polar, double radius) {
		super(smarts);
		this.polar = polar;
		this.radius = radius;
	}
}
