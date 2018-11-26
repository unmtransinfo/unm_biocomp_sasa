/**
 * 
 */
package com.sunset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author oleg
 * clean RDF files
 */
public class RDFClean {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if(args.length != 2) {
			System.out.println("Usage: java RDFClean input.rdf output_clean.rdf");
			return;
		}
		BufferedReader r = new BufferedReader(new FileReader(args[0]));
		PrintWriter w = new PrintWriter(args[1]);
		String line;
		while((line = r.readLine()) != null) {
			if(line.startsWith("$DATUM")) {
				w.println(joinLines(r, line));
			} else {
				w.println(line);
			}
		}
		w.close();
		r.close();
	}
	
	// remove new line characters from $DATUM fields
	
	public static String joinLines(BufferedReader r, String currentLine) throws IOException {
		StringBuilder sb = new StringBuilder(currentLine);
		String line;
		while((line = r.readLine()) != null && !line.startsWith("$")) {
			if(sb.charAt(sb.length() - 1) == '+') {
				sb.append("\n" + line);
			} else {
				sb.append(line);
			}
		}
		if(line != null) {
			return sb.toString() + "\n" + line;
		} else {
			return sb.toString();
		}
	}

}
