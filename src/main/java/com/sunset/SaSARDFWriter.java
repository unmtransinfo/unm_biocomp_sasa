package com.sunset;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import chemaxon.calculations.hydrogenize.Hydrogenize;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.io.MPropHandler;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;

public class SaSARDFWriter {

	private int[] rt = new int[2];
	private int[] pn = new int[2];
	private double[] psa = new double[2];

	private String inputFile = null;
	private String outputFile = null;
	private boolean useVCC = false;
	private boolean verbose = false;
	private String addPropFile = null;

	private VCCLab vccLab = null;

	@SuppressWarnings("static-access")
	public boolean processCommandLine(String[] args) {
		CommandLineParser parser = new PosixParser();
		Options opt = new Options();
		opt.addOption(OptionBuilder.withLongOpt("input-file")
				.isRequired()
				.hasArg()
				.withArgName("Use specified molecular structure input file")
				.create("i"));
		opt.addOption(OptionBuilder.withLongOpt("output-file")
				.isRequired()
				.hasArg()
				.withArgName("Use specified molecular structure output file")
				.create("o"));
		opt.addOption(OptionBuilder.withLongOpt("use-vcc")
				.withArgName("Use VCCLAB to compute ALOGS and ALOGP, default false")
				.create("c"));
		opt.addOption(OptionBuilder.withLongOpt("add-props")
				.isRequired()
				.hasArg()
				.withArgName("Use specified file with additional properties")
				.create("a"));
		opt.addOption(OptionBuilder.withLongOpt("verbose")
				.withArgName("Be verbose")
				.create("v"));
		opt.addOption("h", "help", false, "show this help");
		HelpFormatter helpFormater = new HelpFormatter();
		try {
			CommandLine cmd = parser.parse(opt, args);
			inputFile = cmd.getOptionValue("i");
			outputFile = cmd.getOptionValue("o");
			addPropFile = cmd.getOptionValue("a");
			if(cmd.hasOption("v")) {
				verbose = true;
			}
			if(cmd.hasOption("c")) {
				useVCC = true;
			}

			if(cmd.hasOption("h")) {
				throw new ParseException("Help requested");
			}
		} catch (ParseException ex) {
			System.err.println(ex.getMessage());
			helpFormater.printHelp("SaSARDFWriter", opt);
			return false;
		}
		return true;
	}

	/**
	 * @param args[0] input file
	 * @param args[1] output file
	 */
	public static void main(String[] args) {
		SaSARDFWriter runner = new SaSARDFWriter();
		if(!runner.processCommandLine(args)) {
			return;
		}
		runner.run();
	}

	public void run() {
		MolImporter rdfReader;
		int counter = 1;
		MolExporter rdfWriter;
		BufferedReader addPropsReader;
		Molecule mol = null;
		if(useVCC) {
			vccLab = new VCCLab();
		}
		try {
			rdfReader = new MolImporter(new FileInputStream(inputFile), "rdf");
			rdfWriter = new MolExporter(new FileOutputStream(outputFile), "rdf");
			addPropsReader = new BufferedReader(new FileReader(addPropFile));
			String[] propNames = addPropsReader.readLine().split("\t");
			String[] propValues;
			while((mol = rdfReader.read()) != null) {
				propValues = addPropsReader.readLine().split("\t");
				for(int i = 0; i < propValues.length; i++) {
					if(propValues[i].length() > 0) {
						mol.setProperty(propNames[i], propValues[i]);
					}
				}
				setDescriptors(mol);
				rdfWriter.write(mol);
				if(verbose) {
					System.out.println(counter + " records processed");
				}
				counter++;
			}
			rdfWriter.close();
			rdfReader.close();
			addPropsReader.close();
		} catch (Exception e) {
			System.err.println("ooops error at record: " + counter);
			System.err.println("DB mol record: " + MPropHandler.convertToString(mol.properties(), "ROOT:SMDL.ID"));
			e.printStackTrace();
			return;
		}
	}

	public final void setDescriptors(Molecule mol) throws Exception {
		Molecule copyMol = mol.cloneMolecule();
		copyMol.aromatize(MoleculeGraph.AROM_GENERAL);
		Hydrogenize.removeHAtoms(copyMol);
		mol.setProperty("ROOT:NO.ATOMS", Integer.toString(mol.getAtomCount()));
		mol.setProperty("ROOT:NO.BONDS", Integer.toString(mol.getBondCount()));
		Rings.setMolecule(copyMol);
		int[][] rings = Rings.getSSSR();
		mol.setProperty("ROOT:NO.RINGS", Integer.toString(rings.length));
		RotRigBonds.rotBonds(copyMol, rt);
		mol.setProperty("ROOT:NO.ROT.BONDS", Integer.toString(rt[0]));
		mol.setProperty("ROOT:NO.RIG.BONDS", Integer.toString(mol.getBondCount() - rt[0] - rt[1]));
		AtomCounts.elementCounts(copyMol);
		mol.setProperty("ROOT:NO.HETERO.ATOMS", Integer.toString(AtomCounts.sumHetero()));
		mol.setProperty("ROOT:NO.NONPOL.ATOMS", Integer.toString(AtomCounts.nonPolAtoms()));
		Ionizable.ionizableGroups(copyMol, pn);
		mol.setProperty("ROOT:NO.POS.IONIZ", Integer.toString(pn[0]));
		mol.setProperty("ROOT:NO.NEG.IONIZ", Integer.toString(pn[1]));
		int donnors = HBonds.getDonors(copyMol);
		mol.setProperty("ROOT:LPK.HB.DON", Integer.toString(donnors));
		int acceptors = HBonds.getAcceptors(copyMol);
		mol.setProperty("ROOT:LPK.HB.ACC", Integer.toString(acceptors));
		int lpkScore = 0;
		if(copyMol.getMass() > 500) {
			lpkScore++;
		}
		if(donnors > 5) {
			lpkScore++;
		}
		if(acceptors > 10) {
			lpkScore++;
		}
		double clogp = 0.0d;
		if(MPropHandler.convertToString(mol.properties(), "ROOT:LPK.CLOGP") != null) {
			clogp = Double.parseDouble(MPropHandler.convertToString(mol.properties(), "ROOT:LPK.CLOGP")); 
		} 
		if(clogp > 5) {
			lpkScore++;
		}
		mol.setProperty("ROOT:LPK.SCORE", Integer.toString(lpkScore));
		Hydrogenize.addHAtoms(copyMol);
		PSA.getPSA(copyMol, psa);
		mol.setProperty("ROOT:SFA.POL", Double.toString(psa[0]));
		mol.setProperty("ROOT:SFA.NONPOL", Double.toString(psa[1]));
		mol.setProperty("ROOT:SFA.PERC.POL", Double.toString(psa[0]/(psa[0] + psa[1])));
		mol.setProperty("ROOT:SFA.PERC.NONPOL", Double.toString(psa[1]/(psa[0] + psa[1])));
		Hydrogenize.removeHAtoms(copyMol);
		mol.setProperty("ROOT:MOL.ABE", Double.toString(ABE.calcABE(copyMol)));
		mol.setProperty("ROOT:MOL.SMCM", Double.toString(Complexity.complexity(copyMol)));
		boolean vccSuccess = false;
		if(useVCC && copyMol.getAtomCount() < 255 && AtomCounts.sumCarbons() > 0
				&& AtomCounts.sumMetal() == 0) {
			vccLab.setSMILES(MolExporter.exportToFormat(copyMol, "smiles"));
			try {
				vccLab.run();
				vccSuccess = true;
			} catch (Exception e ) {
				System.err.println("VCC can't calculate properties for: " + MolExporter.exportToFormat(copyMol, "smiles") + "\t" + MPropHandler.convertToString(mol.properties(), "ROOT:SMDL.ID"));
			}
			if(vccSuccess) {
				mol.setProperty("ROOT:EST.LOGKOW", Double.toString(vccLab.getLogP()));
				mol.setProperty("ROOT:EST.LOGWSOL", Double.toString(vccLab.getLogS()));
			}
		}
	}
}
