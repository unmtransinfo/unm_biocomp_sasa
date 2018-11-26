/**
 * 
 */
package com.sunset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Calendar;
import java.util.GregorianCalendar;

/**
 * @author oleg
 * Renumber SMDL.ID, SMDL.IDX, SMDL.SID from a given constant
 *
 */
public class RDFRenumber {

	private static Calendar calendar = new GregorianCalendar();
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if(args.length != 5) {
			System.out.println("Usage: java RDFRenumber RDF_File SDML.ID_start SMDL.SID_start OLD_SMDL.SID RDF_output");
			return;
		}
		// File dir = new File(args[0]);
		
		String[] fileNames = {args[0]};
		int smdl_id = Integer.parseInt(args[1]);
		
		
		int c = smdl_id;
		PrintWriter writer = new PrintWriter(args[4]);
		writer.println("$RDFILE 1\n$DATM " + getDate() + " " + getTime());
		for(String fileName : fileNames) {
			String line;
			BufferedReader reader = new BufferedReader(new FileReader(fileName));
			reader.readLine(); reader.readLine();
			while((line = reader.readLine()) != null) {
				if(line.startsWith("$MFMT $MIREG")) {
					writer.println("$MFMT $MIREG " + c);
					c++;
					continue;
				}

				if(line.equals("$DTYPE ROOT:ID")) {
					writer.println(line);
					reader.readLine();
					writer.println("$DATUM " + smdl_id);
					smdl_id++;
					continue;
				}
				/*if(line.equals("$DTYPE ROOT:SMDL.ID")) {
					writer.println(line);
					reader.readLine();
					writer.println("$DATUM " + smdl_id);
					continue;
				}
				if(line.equals("$DTYPE ROOT:SMDL.IDX")) {
					writer.println(line);
					reader.readLine();
					writer.printf("$DATUM SMDL-%08d%n", smdl_id);
					smdl_id++;
					continue;
				}*/
				/*if(line.equals("$DTYPE ROOT:SMDL.SID")) {
					writer.println(line);
					int s = Integer.parseInt(reader.readLine().split(" ")[1]);
					if(s == ser_sid) {
						writer.println("$DATUM " + sid);
					} else {
						ser_sid = s;
						sid++;
						writer.println("$DATUM " + sid);
					}
					continue;
				}*/
				
				writer.println(line);
			}
			reader.close();
		}
		writer.close();
	}

	public static String getDate() {
		return (calendar.get(Calendar.MONTH) + 1) + "/" + calendar.get(Calendar.DAY_OF_MONTH) + "/" + calendar.get(Calendar.YEAR);
	}
	
	public static String getTime() {
		return calendar.get(Calendar.HOUR) + ":" + calendar.get(Calendar.MINUTE) + ":" + calendar.get(Calendar.SECOND);
	}
	
}
