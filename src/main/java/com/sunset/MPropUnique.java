package com.sunset;

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
 * Filter structures files for unique value of specified properties
 */
public class MPropUnique {
	
	private static class Property {
		String name;
		HashSet<String> values = new HashSet<String>();
		
		public Property(String name) {
			this.name = name;
		}
		
		public boolean filter(Molecule mol) {
			if(mol.properties().get(name) == null) {
				return true;
			}
			String value = MPropHandler.convertToString(mol.properties(), name);
			if(values.contains(value)) {
				return false;
			} else {
				values.add(value);
				return true;
			}
		}	
	}
	
	@Option(name = "-i", aliases = {"--input"}, usage="input structure file", required = true, metaVar = "<FILE>")
	private String inputFileName;
	
	@Option(name = "-o", aliases= {"--output"}, usage="output structure file", required = true, metaVar = "<FILE>")
	private String outputFileName;
	
	@Option(name = "-f", aliases = {"--format"}, required = false, usage = "output file format, default same as input", metaVar = "FORMAT")
	private String fmt;
	
	@Option(name = "-p", aliases = {"--properties"}, required = true, usage = "input property filters", metaVar = "PROP_NAME1,PROP_NAME2,...")
	private List<String> propNames = new ArrayList<String>();
	
	@Option(name = "-h", aliases = {"--help"}, required = false, usage = "show this message")
	private boolean help = false;

	public void run() {
		List<Property> filters = new ArrayList<MPropUnique.Property>();
		for(String propName : propNames) {
			filters.add(new Property(propName));
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
		MPropUnique app = new MPropUnique();
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
