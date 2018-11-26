package com.sunset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import chemaxon.formats.MolFormatException;
import chemaxon.sss.search.MolSearch;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.MoleculeGraph;
import chemaxon.util.MolHandler;


public class SaSA {

	private static String inputFileName = "";
	private static String outFileName = "";
	private static int aromModel = MoleculeGraph.AROM_BASIC;
	
	public static void main(String[] args) throws IOException, MolFormatException, SearchException {
		
		if(!SaSA.parseCommandLine(args)) {
			return;
		}
		
		MolHandler mh = new MolHandler();
		mh.setQueryMode(true);
		
		mh.setMolecule("[!H0;#7,#8;!$([*,-,-2,-3]),!$(*-*=,#*)]");
		MolSearch hbd = new MolSearch();
		hbd.setQuery(mh.getMolecule());
				
		mh.setMolecule("[!$([#6,F,Cl,Br,I,o,$([#8](-[#6a])-[#6a]),s,nX3,$([#7]C=O),$([#7]-[#6a]),#7v5,#15v5,#16,#1,*+1,*+2,*+3])]");
		MolSearch hba = new MolSearch();
		hba.setQuery(mh.getMolecule());
						
		mh.setMolecule("[#8,#7,Pv5,Sv4,Sv6]");
		MolSearch p_at = new MolSearch();
		p_at.setQuery(mh.getMolecule());
						
		mh.setMolecule("[#6,F,Cl,Br,I]");
		MolSearch np_at = new MolSearch();
		np_at.setQuery(mh.getMolecule());
						
		BufferedReader reader = new BufferedReader(new FileReader(inputFileName));
		PrintStream writer = new PrintStream(outFileName);
		writer.println("ID\tP_at\tNP_at\tSHDA\tMWP_at\tMWNP_at\tMWSHDA");
		String line;
		MolHandler mht = new MolHandler();
		int c = 0;
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\\s+");
			writer.print(tokens[1]);
			mht.setMolecule(tokens[0]);
			mht.aromatize(aromModel);
			p_at.setTarget(mht.getMolecule());
			np_at.setTarget(mht.getMolecule());
			hbd.setTarget(mht.getMolecule());
			hba.setTarget(mht.getMolecule());
			int heavy = mht.getHeavyAtomCount();
			double mw = mht.calcMolWeightInDouble();
			int P_at = p_at.getMatchCount();
			int NP_at = np_at.getMatchCount() - P_at;
			int HBA = hba.getMatchCount();
			int HBD = hbd.getMatchCount();
			int SHDA = HBD + HBA;
			double MWP_at = mw * (1D*P_at/heavy);
			double MWNP_at = mw * (1D*NP_at/heavy);
			double MWSHDA = mw * (1D*SHDA/heavy);
			writer.printf("\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n", P_at, NP_at, SHDA, MWP_at, MWNP_at, MWSHDA);
			c++;
			if(c % 1000 == 0) {
				System.out.println(c + " records processed");
			}
		}
		reader.close();
		writer.close();
	}
	
	@SuppressWarnings("static-access")
	public static final boolean parseCommandLine(String[] args) {
		CommandLineParser parser = new PosixParser();
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("input")
				.hasArg()
				.withLongOpt("input")
				.withDescription("use given input SMILES file")
				.isRequired()
				.create("i"));
		
		options.addOption(OptionBuilder.withArgName("output")
				.withLongOpt("output")
				.hasArg()
				.isRequired()
				.withDescription("use given output file")
				.create("o"));
		
		options.addOption(OptionBuilder.withArgName("aromata")
				.withLongOpt("aromata")
				.hasArg()
				.withDescription("aromatization model [b(a)sic default, (g)eneral]")
				.create("a"));
		
		options.addOption("h", "help", false, "show this help");
		
		HelpFormatter helpFormater = new HelpFormatter();
		
		try {
			 
			CommandLine cmd = parser.parse(options, args);
			if(cmd.hasOption("h")) {
				helpFormater.printHelp("SaSA", options);
				return false;
			}
			if(cmd.hasOption("a")) {
				if(cmd.getOptionValue("a").equals("g")) {
					aromModel = MoleculeGraph.AROM_GENERAL;
				} else if (cmd.getOptionValue("a").equals("a")) {
					aromModel = MoleculeGraph.AROM_BASIC;
				} else {
					helpFormater.printHelp("SaSA", options);
					return false;
				}
				
			}
			inputFileName = cmd.getOptionValue("i");
			outFileName = cmd.getOptionValue("o");
		} catch (ParseException ex) {
			System.out.println(ex.getMessage());
			helpFormater.printHelp("SaSA", options);
			return false;
		}
		if(inputFileName == null || outFileName == null) {
			helpFormater.printHelp("SaSA", options);
			return false;
		}
		return true;
	}

	
}
