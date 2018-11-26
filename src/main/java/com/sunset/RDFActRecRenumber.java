package com.sunset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

/**
 * @author Oleg Ursu
 * Renumber activity records
 */
public class RDFActRecRenumber {

	private static final Logger logger = LoggerFactory.getLogger(RDFActRecRenumber.class);
	
	private boolean help;
	
	@Parameter(names = {"-i", "--input"}, required = true, description = "input RDF file")
	private String inputRDFFileName;
	
	@Parameter(names = {"-o", "--output"}, required = true, description = "output RDF file")
	private String outputRDFFileName;
	
	public boolean help() {
		return this.help;
	}
	
	public void run() {
		BufferedReader reader = null;
		PrintWriter writer = null;
		try {
			reader = new BufferedReader(new FileReader(inputRDFFileName));
			writer = new PrintWriter(outputRDFFileName);
			String line;
			int ridx = 0;
			long lineCount = 0L;
			int wrongidx = 0;
			boolean wrongBlock = false;
			while((line = reader.readLine()) != null) {
				lineCount++;
				if(line.startsWith("$MFMT $MIREG ")) {
					ridx = 0;
					wrongidx = 0;
					writer.println(line);
					wrongBlock = false;
					continue;
				}
				if(line.startsWith("$DTYPE ROOT:ACT.LIST(")) {
					String str = line.substring(21, line.indexOf(')'));
					int i = Integer.parseInt(str);
					if(!wrongBlock && (i == ridx || i == (ridx + 1))) {
						ridx = i;
						writer.println(line);
					} else {
						wrongBlock = true;
						if(i == wrongidx) {
							writer.println("$DTYPE ROOT:ACT.LIST(" + ridx + line.substring(line.indexOf(')')));
						} else {
							ridx++;
							wrongidx = i;
							writer.println("$DTYPE ROOT:ACT.LIST(" + ridx + line.substring(line.indexOf(')')));
							// replace $DATUM XXX line
							reader.readLine();
							writer.println("$DATUM " + ridx);
						}
						logger.info("detected misnumber at {}", lineCount);
						logger.info("line string: {}", line);
					}
				} else {
					writer.println(line);
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if(writer != null) {
				writer.close();
			}
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		RDFActRecRenumber app = new RDFActRecRenumber();
		JCommander jcmdr = new JCommander(app, args);
		if(app.help) {
			jcmdr.usage();
			return;
		}
		app.run();
	}

}
