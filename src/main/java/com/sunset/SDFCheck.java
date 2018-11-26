package com.sunset;
import java.io.FileInputStream;
import java.io.FileOutputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.license.LicenseException;
import chemaxon.reaction.Standardizer;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;


public class SDFCheck {

	private static Standardizer molstd;
	private String inputFileName;
	private String outputFileName;
	private String failedFileName;
	private boolean verbose = false;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		SDFCheck checker = new SDFCheck("keepone..neutralize..dearomatize..dehydrogenize..O=N=O>>[O-]-[N+]=O..N=[N:1]#[N:2]>>N=[N+:1]=[N-:2]");
		if(!checker.processCommandLine(args)) {
			return;
		}
		checker.run();
	}
	
	public SDFCheck(String stdConfig) throws Exception {
		if(molstd == null) {
			molstd = new Standardizer(stdConfig);
		}
	}
	
	public void standardize(Molecule mol) throws SearchException, LicenseException {
		molstd.standardize(mol);
	}
	
	public boolean checkCharge(Molecule mol) {
		int charge = mol.getFormalCharge();
		if(charge != 0) {
			for(MolAtom atom : mol.getAtomArray()) {
				if(atom.getCharge() == charge && atom.getAtno() != 7) {
					return false;
				}
			}
		}
		return true;
	}
	
	public void run() throws Exception {
		MolImporter mreader = new MolImporter(new FileInputStream(inputFileName), "sdf");
		MolExporter mwriter = new MolExporter(new FileOutputStream(outputFileName), "sdf");
		MolExporter failedMolWriter = new MolExporter(new FileOutputStream(failedFileName), "sdf");
		Molecule mol;
		int count = 0;
		while((mol = mreader.read()) != null) {
			standardize(mol);
			if(mol.getFragCount(MoleculeGraph.FRAG_KEEPING_MULTICENTERS) > 1) {
				failedMolWriter.write(mol);
				continue;
			}
			if(!checkCharge(mol)) {
				failedMolWriter.write(mol);
				continue;
			}
			if(mol.hasValenceError()) {
				failedMolWriter.write(mol);
				continue;
			}
			if(mol.isQuery()) {
				failedMolWriter.write(mol);
				continue;
			}
			mwriter.write(mol);
			count++;
			if(verbose && count % 100 == 0) {
				System.out.println(count + " records processed");
			}
		}
		mwriter.close();
		failedMolWriter.close();
		mreader.close();
	}
	
	@SuppressWarnings("static-access")
	public boolean processCommandLine(String[] args) {
		CommandLineParser parser = new PosixParser();
		Options opt = new Options();
		opt.addOption(OptionBuilder.withLongOpt("input-file")
									.isRequired()
									.hasArg()
									.withArgName("Use specified SDF input file")
									.create("i"));
		opt.addOption(OptionBuilder.withLongOpt("output-file")
									.isRequired()
									.hasArg()
									.withArgName("Use specified SDF output file")
									.create("o"));
		opt.addOption(OptionBuilder.withLongOpt("bad-file")
									.isRequired()
									.hasArg()
									.withArgName("Use specified SDF output file to write bad structures")
									.create("b"));
		opt.addOption(OptionBuilder.withLongOpt("verbose")
				   					.withArgName("Be verbose")
				   					.create("v"));
		opt.addOption("h", "help", false, "show this help");
		HelpFormatter helpFormater = new HelpFormatter();
		try {
			CommandLine cmd = parser.parse(opt, args);
			if(cmd.hasOption("i")) {
				inputFileName = cmd.getOptionValue("i");
			}
			if(cmd.hasOption("o")) {
				outputFileName = cmd.getOptionValue("o");
			}
			if(cmd.hasOption("b")) {
				failedFileName = cmd.getOptionValue("b");
			}
			if(cmd.hasOption("v")) {
				verbose = true;
			}
			if(cmd.hasOption("h")) {
				throw new ParseException("Help requested");
			}
		} catch (ParseException ex) {
			System.err.println(ex.getMessage());
			helpFormater.printHelp("SDFCheck", opt);
			return false;
		}
		return true;
	}

}
