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
 * Renumber RDF reference file
 */
public class RDFRefRenumber {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if(args.length != 3) {
			System.out.println("Usage: java RDFRefRenumber RDF_Ref_input SMDL_SID_start RDF_Ref_output");
			return;
		}
		BufferedReader reader = new BufferedReader(new FileReader(args[0]));
		int sid = Integer.parseInt(args[1]);
		PrintWriter writer = new PrintWriter(args[2]);
		writer.println(reader.readLine());
		writer.println(reader.readLine());
		String line;
		while((line = reader.readLine()) != null) {
			if(line.startsWith("$RIREG")) {
				writer.println("$RIREG " + sid);
				continue;
			}
			if(line.equals("$DTYPE ROOT:ID")) {
				writer.println(line);
				reader.readLine();
				writer.println("$DATUM " + sid);
				writer.println("$DTYPE ROOT:SMDL_SID");
				writer.println("$DATUM " + sid);
				sid++;
				continue;
			}
			writer.println(line);
		}
		reader.close();
		writer.close();
	}

}
