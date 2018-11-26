package com.sunset;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.io.MPropHandler;
import chemaxon.struc.Molecule;

public class RDFAddIds {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if(args.length != 4) {
			System.out.println("Usage: java RDFAddIds rdf_input ref_ids rdf_out start");
			return;
		}
		String inputFileName = args[0];
		String refInputFileName = args[1];
		String outputFileName = args[2];
		int id = Integer.parseInt(args[3]);
		HashMap<String, Integer> refs = new HashMap<String, Integer>();
		BufferedReader refReader = new BufferedReader(new FileReader(refInputFileName));
		String line;
		String[] tokens;
		while((line = refReader.readLine()) != null) {
			tokens = line.split("\t");
			refs.put(tokens[1], new Integer(tokens[0]));
		}
		refReader.close();
		MolImporter rdfReader = new MolImporter(new FileInputStream(inputFileName), "rdf");
		MolExporter rdfWriter = new MolExporter(new FileOutputStream(outputFileName), "rdf");
		Molecule mol;
		String ref;
		int smdl_sid;
		while((mol = rdfReader.read()) != null) {
			mol.setProperty("ROOT:ID", Integer.toString(id));
			mol.setProperty("ROOT:SMDL.ID", Integer.toString(id));
			mol.setProperty("ROOT:SMDL.IDX", String.format("SMDL-%08d", id));
			ref = MPropHandler.convertToString(mol.properties(), "ROOT:MOL.REF");
			smdl_sid = refs.get(ref);
			mol.setProperty("ROOT:SMDL.SID", Integer.toString(smdl_sid));
			rdfWriter.write(mol);
			id++;
		}
		rdfReader.close();
		rdfWriter.close();
	}
}
