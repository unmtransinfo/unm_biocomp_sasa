/**
 * 
 */
package com.sunset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Concatenate RDF file lines in $DATUM fields longer than 80 characters.
 * @author Oleg Ursu
 *
 */
public class RDFConcatLines {

	/**
	 * @param args input file in args[0]
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if(args.length != 1) {
			System.out.println("Usage: java RDFConcatLines RDF_file");
			return;
		}
		BufferedReader reader = new BufferedReader(new FileReader(args[0]));
		String line;
		int l;
		while((line = reader.readLine()) != null) {
			if(line.length() == 80 + 1 && line.charAt(80) == '+') {
				l = 0;
				while(line.length()==(++l*80 + 1) && '+'== line.charAt(l*80)){
	                line=line.substring(0,l*80)+reader.readLine();
	            }
				System.out.println(line);
			} else {
				System.out.println(line);
			}
		}
		System.out.flush();
		reader.close();
	}

}
