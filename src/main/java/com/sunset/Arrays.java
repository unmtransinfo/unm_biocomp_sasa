/**
 * 
 */
package com.sunset;

/**
 * @author oleg
 *
 */
public class Arrays {
	public static int indexOf(int i, int[] a) {
		for(int j = 0; j < a.length; j++) {
			if(a[j] == i) {
				return j;
			}
		}
		return -1;
	}
}
