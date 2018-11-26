/**
 * 
 */
package com.sunset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import chemaxon.descriptors.Metrics;
import chemaxon.formats.MolImporter;
import chemaxon.sss.search.MolSearch;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;

/**
 * @author Oleg Ursu
 * compare structures in two input SMILES files based on SMDL descriptors
 * output 1 if tanimoto = 1, 0 otherwise
 */
public class SMDLCompare {

	private static final ArrayList<MolSearch> SMDL = new ArrayList<MolSearch>();
	private static final ArrayList<Molecule> QUERY = new ArrayList<Molecule>();
	private static final ArrayList<String> QNAME = new ArrayList<String>();
	private static final ArrayList<int[]> QSMDL = new ArrayList<int[]>();
	private static final ArrayList<String> QSMILES = new ArrayList<String>();
	
	
	/**
	 * @param args[0] = SMDL descriptors file, args[1] = query file, args[2] = target file
	 */
	public static void main(String[] args) throws IOException, SearchException {
		if(args.length != 3) {
			System.out.println("usage: java SMDLCompare descriptors query target");
			return;
		}
		loadSMDL(args[0]);
		
		loadQuery(args[1]);
		
		String line;
		String[] tokens;
		BufferedReader reader = new BufferedReader(new FileReader(args[2]));
		boolean found;
		int[] smdl;
		Molecule mol;
		System.out.println("mol.smi\tID\tmol.isdrug_new\tdrug_name");
		while((line = reader.readLine()) != null) {
			tokens = line.split("\t");
			found = false;
			for(int i = 0; i < QSMILES.size(); i++) {
				if(tokens[0].equals(QSMILES.get(i))) {
					System.out.println(line + "\t" + 1 + "\t" + QNAME.get(i));
					found = true;
					break;
				}
			}
			if(!found && tokens.length >= 4) {
				for(int i = 0; i < QNAME.size(); i++) {
					if(QNAME.get(i).compareToIgnoreCase(cleanName(tokens[3])) == 0) {
						System.out.println(line + "\t" + 1 + "\t" + QNAME.get(i));
						found = true;
						break;
					}
				}
				
			}
			if(!found) {
				mol = MolImporter.importMol(tokens[0], "smiles");
				mol.aromatize(MoleculeGraph.AROM_GENERAL);
				smdl = getSMDL(mol);
				for(int i = 0; i < QSMDL.size(); i++) {
					if(Metrics.tanimoto(smdl, QSMDL.get(i)) == 0) {
						System.out.println(line + "\t" + 1 + "\t" + QNAME.get(i));
						found = true;
						break;
					}
				}
			}
			if(!found) {
				System.out.println(line + "\t" + 0 + "\t");
			}
		}
		reader.close();
	}
	
	public static final void loadSMDL(String fileName) throws IOException {
		SMDL.clear();
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		String line;
		while((line = reader.readLine()) != null) {
			if(line.startsWith("#")) {
				continue;
			}
			MolSearch ms = new MolSearch();
			ms.setQuery(MolImporter.importMol(line, "smarts:d"));
			SMDL.add(ms);
		}
		reader.close();
	}
	
	public static final void loadQuery(String fileName) throws IOException, SearchException {
		QUERY.clear();
		QNAME.clear();
		QSMDL.clear();
		QSMILES.clear();
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		String line;
		String[] tokens;
		Molecule mol;
		while((line = reader.readLine()) != null) {
			tokens = line.split("\t");
			mol = MolImporter.importMol(tokens[0], "smiles");
			mol.aromatize(MoleculeGraph.AROM_GENERAL);
			QUERY.add(mol);
			QNAME.add(cleanName(tokens[1]));
			QSMDL.add(getSMDL(mol));
			QSMILES.add(tokens[0]);
		}
		reader.close();
	}
	
	public static final int[] getSMDL(Molecule mol) throws SearchException {
		int[] smdl = new int[SMDL.size()];
		MolSearch ms;
		int[][] matches;
		for(int i = 0; i < SMDL.size(); i++) {
			ms = SMDL.get(i);
			ms.setTarget(mol);
			matches = ms.findAll();
			if(matches != null) {
				smdl[i] = matches.length;
			} else {
				smdl[i] = 0;
			}
		}
		return smdl;
	}
	
	public static final String cleanName(String name) {
		String cname = name.trim();
		cname = cname.replaceAll("\"", "");
		return cname;
	}
	
}
