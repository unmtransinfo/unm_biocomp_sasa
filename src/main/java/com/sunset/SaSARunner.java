package com.sunset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

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
import chemaxon.marvin.calculations.logPPlugin;
import chemaxon.marvin.io.MPropHandler;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;

/**
 * @author Oleg Ursu
 * simple command line interface to compute SaSA descriptors from a SMILES file
 *
 */
public class SaSARunner {

	private int[] rt = new int[2];
	private int[] pn = new int[2];
	private double[] psa = new double[2];

	private String inputFile = null;
	private String outputFile = null;
	private boolean useVCC = false;
	private String sdTagName = null;
	private boolean verbose = false;
	private String addPropFileName = null;
	private String delim = "\t";

	private String fmt = "smiles:T*";

	private VCCLab vccLab = null;

	private logPPlugin logp;

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
				.withArgName("Use specified file to write descriptors")
				.create("o"));
		opt.addOption(OptionBuilder.withLongOpt("identifier")
				.hasArg()
				.withArgName("Use the following sourse for molecule identifier, default molecule title if present, else sequental integer value")
				.create("e"));
		opt.addOption(OptionBuilder.withLongOpt("format")
				.hasArg()
				.withArgName("Use specified structure format in output file, default SMILES")
				.create("f"));
		opt.addOption(OptionBuilder.withLongOpt("use_vcc")
				.withArgName("Use VCCLAB to compute ALOGS and ALOGP, default false")
				.create("c"));
		opt.addOption(OptionBuilder.withLongOpt("add-propperties")
				.hasArg()
				.withArgName("Read additional properties from file, default NONE")
				.create("a"));
		opt.addOption(OptionBuilder.withLongOpt("propperties-delim")
				.hasArg()
				.withArgName("Properties file delimiter, default TAB")
				.create("d"));
		opt.addOption(OptionBuilder.withLongOpt("verbose")
				.withArgName("Be verbose")
				.create("v"));
		opt.addOption("h", "help", false, "show this help");
		HelpFormatter helpFormater = new HelpFormatter();
		try {
			CommandLine cmd = parser.parse(opt, args);
			if(cmd.hasOption("i")) {
				inputFile = cmd.getOptionValue("i");
			}
			if(cmd.hasOption("o")) {
				outputFile = cmd.getOptionValue("o");
			}
			if(cmd.hasOption("e")) {
				sdTagName = cmd.getOptionValue("e");
			}
			if(cmd.hasOption("v")) {
				verbose = true;
			}
			if(cmd.hasOption("c")) {
				useVCC = true;
			}
			if(cmd.hasOption("f")) {
				fmt = cmd.getOptionValue("f");
			}
			if(cmd.hasOption("a")) {
				addPropFileName = cmd.getOptionValue("a");
			}
			if(cmd.hasOption("d")) {
				delim = cmd.getOptionValue("d");
			}
			if(cmd.hasOption("h")) {
				throw new ParseException("Help requested");
			}
		} catch (ParseException ex) {
			System.err.println(ex.getMessage());
			helpFormater.printHelp("SaSARunner", opt);
			return false;
		}
		return true;
	}

	/**
	 * @param args[0] input file
	 * @param args[1] output file
	 */
	public static void main(String[] args) {
		SaSARunner runner = new SaSARunner();
		if(!runner.processCommandLine(args)) {
			return;
		}
		runner.run();
	}

	public void run() {
		MolImporter reader;
		int counter = 1;
		MolExporter writer;
		Molecule mol;
		String molId;
		if(useVCC) {
			vccLab = new VCCLab();
		}
		BufferedReader addPropReader = null;
		ArrayList<String> propNames = null;
		try {
			logp = new logPPlugin();
			logp.setlogPMethod(logPPlugin.METHOD_WEIGHTED);
			logp.setWeightOfMethods(1, 2, 1, 0);
			logp.setCloridIonConcentration(0.2);
			logp.setNaKIonConcentration(0.2);
			logp.setUserTypes("logPTrue");
			String line;
			if(addPropFileName != null) {
				addPropReader = new BufferedReader(new FileReader(addPropFileName));
				line = addPropReader.readLine();
				propNames = new ArrayList<String>();
				for(String propName : line.split(delim)) {
					propNames.add(propName);
				}
			}
			reader = new MolImporter(inputFile);
			writer = new MolExporter(outputFile, fmt);
			while((mol = reader.read()) != null) {
				if(sdTagName != null) {
					molId = MPropHandler.convertToString(mol.properties(), sdTagName);
				} else if(mol.getName().length() > 0) {
					molId = mol.getName();
				} else {
					molId = Integer.toString(counter);
				}
				mol.setName(molId);
				calculateDescriptors(mol);
				if(addPropReader != null) {
					line = addPropReader.readLine();
					String[] tokens = line.split(delim);
					for(int i = 0; i < tokens.length; i++) {
						if(tokens[i].length() > 0) {
							mol.setProperty(propNames.get(i), tokens[i]);
						}
					}
				}
				logp.setMolecule(mol);
				logp.run();
				double logPValue = logp.getlogPTrue();
				if(Double.compare(logPValue, Double.NaN) != 0) {
					mol.setProperty("logP", Double.toString(logPValue));
				}
				writer.write(mol);
				if(verbose) {
					System.out.println(counter + " records processed");
				}
				counter++;
			}
			writer.close();
			reader.close();
			if(addPropReader != null) {
				addPropReader.close();
			}
		} catch (Exception e) {
			System.err.println("ooops error at record: " + counter);
			return;
		}
	}

	/*	public final String getDescriptors(Molecule mol) throws Exception {
		StringBuilder sb = new StringBuilder();
		Hydrogenize.removeHAtoms(mol);
		sb.append("\t" + mol.getMass() + delimiter + mol.getAtomCount() + delimiter + mol.getBondCount());
		Rings.setMolecule(mol);
		int[][] rings = Rings.getSSSR();
		sb.append(delimiter + rings.length);
		RotRigBonds.rotBonds(mol, rt);
		sb.append(delimiter + rt[0] + delimiter + (mol.getBondCount() - rt[0] - rt[1]));
		AtomCounts.elementCounts(mol);
		sb.append(delimiter + AtomCounts.sumHetero() + delimiter + AtomCounts.nonPolAtoms());
		Ionizable.ionizableGroups(mol, pn);
		sb.append(delimiter + pn[0] + delimiter + pn[1]);
		sb.append(delimiter + HBonds.getDonors(mol) + delimiter + HBonds.getAcceptors(mol));
		Hydrogenize.addHAtoms(mol);
		PSA.getPSA(mol, psa);
		sb.append(delimiter + psa[0] + delimiter + psa[1] + delimiter + (psa[0]/(psa[0] + psa[1])) + delimiter + (psa[1]/(psa[0] + psa[1])));
		Hydrogenize.removeHAtoms(mol);
		sb.append(delimiter + ABE.calcABE(mol));
		sb.append(delimiter + Complexity.complexity(mol));
		sb.append(delimiter + (Rings.aliphaticRings() != null ? Rings.aliphaticRings().length : "0"));
		sb.append(delimiter + (Rings.aromaticRings() != null ? Rings.aromaticRings().length : "0"));
		sb.append(delimiter + (Rings.heteroAliphaticRings() != null ? Rings.heteroAliphaticRings().length : "0"));
		sb.append(delimiter + (Rings.heteroAromaticRings() != null ? Rings.heteroAromaticRings().length : "0"));
		boolean vccSuccess = false;
		if(useVCC && mol.getAtomCount() < 255 && AtomCounts.sumCarbons() > 0) {
			vccLab.setSMILES(MolExporter.exportToFormat(mol, "smiles"));
			try {
				vccLab.run();
				vccSuccess = true;
			} catch (Exception e ) {
				System.err.println("VCC can't calculate properties for: " + MolExporter.exportToFormat(mol, "smiles"));
			}
			if(vccSuccess) {
				sb.append(delimiter + vccLab.getLogP() + delimiter + vccLab.getLogS());
			}
		}
		return sb.toString();
	}
	 */	
	public final void calculateDescriptors(Molecule smol) throws Exception {
		Molecule mol = smol.cloneMolecule();
		mol.aromatize(MoleculeGraph.AROM_GENERAL);
		Hydrogenize.removeHAtoms(mol);
		smol.setProperty("MOL_WEIGHT", Double.toString(mol.getMass()));
		smol.setProperty("NO_ATOMS", Integer.toString(mol.getAtomCount()));
		smol.setProperty("NO_BONDS", Integer.toString(mol.getBondCount()));
		Rings.setMolecule(mol);
		int[][] rings = Rings.getSSSR();
		smol.setProperty("NO_RINGS", Integer.toString(rings.length));
		RotRigBonds.rotBonds(mol, rt);
		smol.setProperty("NO_ROT_BONDS", Integer.toString(rt[0]));
		smol.setProperty("NO_RIG_BONDS", Integer.toString(mol.getBondCount() - rt[0] - rt[1]));
		AtomCounts.elementCounts(mol);
		smol.setProperty("NO_HETERO_ATOMS", Integer.toString(AtomCounts.sumHetero()));
		smol.setProperty("NO_NONPOL_ATOMS", Integer.toString(AtomCounts.nonPolAtoms()));
		Ionizable.ionizableGroups(mol, pn);
		smol.setProperty("NO_POS_IONIZ", Integer.toString(pn[0]));
		smol.setProperty("NO_NEG_IONIZ", Integer.toString(pn[1]));
		smol.setProperty("LPK_HB_DON", Integer.toString(HBonds.getDonors(mol)));
		smol.setProperty("LPK_HB_ACC", Integer.toString(HBonds.getAcceptors(mol)));
		Hydrogenize.addHAtoms(mol);
		PSA.getPSA(mol, psa);
		smol.setProperty("SFA_POL", Double.toString(psa[0]));
		smol.setProperty("SFA_NONPOL", Double.toString(psa[1]));
		smol.setProperty("SFA_PERC_POL", Double.toString(psa[0]/(psa[0] + psa[1])));
		smol.setProperty("SFA_PERC_NONPOL", Double.toString(psa[1]/(psa[0] + psa[1])));
		smol.setProperty("VDW_VOLUME", Double.toString(VABC.calculate(mol)));
		Hydrogenize.removeHAtoms(mol);
		smol.setProperty("MOL_ABE", Double.toString(ABE.calcABE(mol)));
		smol.setProperty("MOL_SMCM", Double.toString(Complexity.complexity(mol)));
		smol.setProperty("NO_ALIPHATIC_RINGS", Integer.toString((Rings.aliphaticRings() != null ? Rings.aliphaticRings().length : 0)));
		smol.setProperty("NO_AROMATIC_RINGS", Integer.toString((Rings.aromaticRings() != null ? Rings.aromaticRings().length : 0)));
		smol.setProperty("NO_HETEROALIPHATIC_RINGS", Integer.toString((Rings.heteroAliphaticRings() != null ? Rings.heteroAliphaticRings().length : 0)));
		smol.setProperty("NO_HETEROAROMATIC_RINGS", Integer.toString((Rings.heteroAromaticRings() != null ? Rings.heteroAromaticRings().length : 0)));
		boolean vccSuccess = false;
		if(useVCC && mol.getAtomCount() < 255 && AtomCounts.sumCarbons() > 0) {
			vccLab.setSMILES(MolExporter.exportToFormat(mol, "smiles"));
			try {
				vccLab.run();
				vccSuccess = true;
			} catch (Exception e ) {
				System.err.println("VCC can't calculate properties for: " + MolExporter.exportToFormat(mol, "smiles"));
			}
			if(vccSuccess) {
				smol.setProperty("ALOGP", Double.toString(vccLab.getLogP()));
				smol.setProperty("ALOGS", Double.toString(vccLab.getLogS()));
			}
		}
	}
}
