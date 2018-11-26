/**
 * 
 */
package com.sunset;

/**
 * @author Oleg Ursu
 *
 */
public class HBPattern extends Pattern {
	public int factor;
	public HBPattern(String smarts, int factor) {
		super(smarts);
		this.factor = factor;
	}
}
