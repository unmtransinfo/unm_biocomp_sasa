package com.sunset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.kohsuke.args4j.ParserProperties;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.io.MPropHandler;
import chemaxon.struc.Molecule;

/**
 * @author Oleg Ursu
 * 
 * Filter structure files using specified property filters
 * 
 */

public class MPropFilter {
	
	private static class Property {
		String name;
		boolean inverse = false;
		HashSet<String> values = new HashSet<String>();
		
		public Property(String fileName) {
			BufferedReader reader;
			try {
				reader = new BufferedReader(new FileReader(fileName));
				String line = reader.readLine();
				String[] tokens = line.split("\\s+");
				this.name = tokens[0];
				this.inverse = Boolean.parseBoolean(tokens[1]);
				while((line = reader.readLine()) != null) {
					values.add(line);
				}
				reader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		public boolean filter(Molecule mol) {
			if(mol.properties().get(name) == null) {
				return true;
			}
			String value = MPropHandler.convertToString(mol.properties(), name);
			if(values.contains(value)) {
				if(inverse) {
					return false;
				} else {
					return true;
				}
			} else {
				if(inverse) {
					return true;
				} else {
					return false;
				}
			}
		}		
	}
	
	@Option(name = "-i", aliases = {"--input"}, usage="input structure file", required = true, metaVar = "<FILE>")
	private String inputFileName;
	
	@Option(name = "-o", aliases= {"--output"}, usage="output structure file", required = true, metaVar = "<FILE>")
	private String outputFileName;
	
	@Option(name = "-f", aliases = {"--format"}, required = false, usage = "output file format, default same as input", metaVar = "FORMAT")
	private String fmt;
	
	@Option(name = "-p", aliases = {"--properties"}, required = true, usage = "input property filters", metaVar = "PROP_FILE1,PROP_FILE2,...")
	private List<String> propFileNames = new ArrayList<String>();
	
	@Option(name = "-h", aliases = {"--help"}, required = false, usage = "show this message")
	private boolean help = false;
	
	public void run() {
		List<Property> filters = new ArrayList<MPropFilter.Property>();
		for(String fileName : propFileNames) {
			filters.add(new Property(fileName));
		}
		MolImporter reader;
		try {
			reader = new MolImporter(inputFileName);
			if(fmt == null) {
				fmt = reader.getFormat();
			}
			MolExporter writer = new MolExporter(outputFileName, fmt);
			Molecule mol;
			boolean pass;
			while((mol = reader.read()) != null) {
				pass = true;
				for(Property filter : filters) {
					if(!filter.filter(mol)) {
						pass = false;
						break;
					}
				}
				if(pass) {
					writer.write(mol);
				}
			}
			writer.close();
			reader.close();
		} catch (MolFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		MPropFilter app = new MPropFilter();
		ParserProperties properties = ParserProperties.defaults();
        properties.withUsageWidth(200);
        CmdLineParser parser = new CmdLineParser(app, properties);
		try {
			parser.parseArgument(args);
			if(app.help()) {
				parser.printUsage(System.out);
				return;
			}
			app.run();
		} catch(CmdLineException ex) {
			parser.printUsage(System.out);
		}
	}

	public boolean help() {
		return help;
	}
}
