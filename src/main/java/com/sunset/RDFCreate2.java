/**
 * 
 */
package com.sunset;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.io.MPropHandler;
import chemaxon.struc.Molecule;

/**
 * @author Oleg Ursu
 * create an hierarchical RDF file from SDF and tab files
 */
public class RDFCreate2 {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if(args.length != 1) {
			System.err.println("Usage: java RDFCreate2 output.rdf");
			return;
		}
		MolImporter molReader = new MolImporter(new FileInputStream("STRUCTURES.SDF"), "sdf");
		Molecule mol;
		HashMap<String, Molecule> molMap = new HashMap<String, Molecule>();
		HashMap<String, Integer> countMap = new HashMap<String, Integer>();
		while((mol = molReader.read()) != null) {
			Molecule cloneMol = mol.cloneMolecule();
			cloneMol.properties().clear();
			for(String propName : mol.properties().getKeys()) {
				cloneMol.setProperty("ROOT:" + propName, MPropHandler.convertToString(mol.properties(), propName));
			}
			cloneMol.setName(MPropHandler.convertToString(mol.properties(), "MOLNAME"));
			molMap.put(MPropHandler.convertToString(mol.properties(), "ID"), cloneMol);
			countMap.put(MPropHandler.convertToString(mol.properties(), "ID"), Integer.valueOf(1));
		}
		molReader.close();
		BufferedReader reader = new BufferedReader(new FileReader("DRUG_TARGET.tab"));
		String[] tags = reader.readLine().split("\t");
		String line;
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			mol = molMap.get(tokens[0]); 
			for(int i = 1; i < tokens.length; i++) {
				String propValue = tokens[i].replace("\"", "");
				if(propValue.length() > 0) {
					mol.setProperty("ROOT:DRUG_TARGET(" + countMap.get(tokens[0]) + "):" + tags[i], propValue);
				}
			}
			countMap.put(tokens[0], Integer.valueOf(countMap.get(tokens[0]) + 1));
		}
		reader.close();
		resetCountMap(countMap);
		
		reader = new BufferedReader(new FileReader("HERG_BIOACT.tab"));
		tags = reader.readLine().split("\t");
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			mol = molMap.get(tokens[0]);
			for(int i = 1; i < tokens.length; i++) {
				String propValue = tokens[i].replace("\"", "");
				if(propValue.length() > 0) {
					mol.setProperty("ROOT:HERG_BIOACT(" + countMap.get(tokens[0]) + "):" + tags[i], propValue);
				}
			}
			countMap.put(tokens[0], Integer.valueOf(countMap.get(tokens[0]) + 1));
		}
		reader.close();
		resetCountMap(countMap);
		
		HashMap<String, String> mixt2StructMap = new HashMap<String, String>();
		HashMap<String, Integer> mixtSeqMap = new HashMap<String, Integer>();
		HashMap<String, Integer> mixtCountMap = new HashMap<String, Integer>();
		
		reader = new BufferedReader(new FileReader("MIXTURE.tab"));
		tags = reader.readLine().split("\t");
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			mol = molMap.get(tokens[0]);
			mixt2StructMap.put(tokens[1], tokens[0]);
			mixtCountMap.put(tokens[1], Integer.valueOf(1));
			for(int i = 2; i < tokens.length; i++) {
				String propValue = tokens[i].replace("\"", "");
				if(propValue.length() > 0) {
					mol.setProperty("ROOT:MIXTURE(" + countMap.get(tokens[0]) + "):" + tags[i], propValue);
				}
			}
			mixtSeqMap.put(tokens[1], countMap.get(tokens[0]));
			countMap.put(tokens[0], Integer.valueOf(countMap.get(tokens[0]) + 1));
		}
		reader.close();
		resetCountMap(countMap);
		
		reader = new BufferedReader(new FileReader("MIXTURE_LD50.tab"));
		tags = reader.readLine().split("\t");
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			mol = molMap.get(mixt2StructMap.get(tokens[0]));
			Integer mixtSeq = mixtSeqMap.get(tokens[0]);
			for(int i = 1; i < tokens.length; i++) {
				String propValue = tokens[i].replace("\"", "");
				if(propValue.length() > 0) {
					mol.setProperty("ROOT:MIXTURE(" + mixtSeq + "):LD50(" + mixtCountMap.get(tokens[0]) + "):" + tags[i], propValue);
				}
			}
			mixtCountMap.put(tokens[0], Integer.valueOf(mixtCountMap.get(tokens[0]) + 1));
		}
		reader.close();
		resetCountMap(mixtCountMap);
		
		reader = new BufferedReader(new FileReader("MIXTURE_WATER_SOL.tab"));
		tags = reader.readLine().split("\t");
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			mol = molMap.get(mixt2StructMap.get(tokens[0]));
			Integer mixtSeq = mixtSeqMap.get(tokens[0]);
			for(int i = 1; i < tokens.length; i++) {
				String propValue = tokens[i].replace("\"", "");
				if(propValue.length() > 0) {
					mol.setProperty("ROOT:MIXTURE(" + mixtSeq + "):WATER_SOLUBILITY(" + mixtCountMap.get(tokens[0]) + "):" + tags[i], propValue);
				}
			}
			mixtCountMap.put(tokens[0], Integer.valueOf(mixtCountMap.get(tokens[0]) + 1));
		}
		reader.close();
		resetCountMap(mixtCountMap);
		
		reader = new BufferedReader(new FileReader("PHASE1_METAB_ENZYME.tab"));
		tags = reader.readLine().split("\t");
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			mol = molMap.get(tokens[0]);
			for(int i = 1; i < tokens.length; i++) {
				String propValue = tokens[i].replace("\"", "");
				if(propValue.length() > 0) {
					mol.setProperty("ROOT:PHASE1_METAB_ENZYME(" + countMap.get(tokens[0]) + "):" + tags[i], propValue);
				}
			}
			countMap.put(tokens[0], Integer.valueOf(countMap.get(tokens[0]) + 1));
		}
		reader.close();
		resetCountMap(countMap);
		
		reader = new BufferedReader(new FileReader("WOMBAT_ACTLIST.tab"));
		tags = reader.readLine().split("\t");
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			mol = molMap.get(tokens[0]);
			for(int i = 1; i < tokens.length; i++) {
				String propValue = tokens[i].replace("\"", "");
				if(propValue.length() > 0) {
					mol.setProperty("ROOT:WOMBAT(" + countMap.get(tokens[0]) + "):" + tags[i], propValue);
				}
			}
			countMap.put(tokens[0], Integer.valueOf(countMap.get(tokens[0]) + 1));
		}
		reader.close();
		resetCountMap(countMap);
		
		reader = new BufferedReader(new FileReader("SALES.tab"));
		tags = reader.readLine().split("\t");
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			mol = molMap.get(tokens[0]);
			for(int i = 1; i < tokens.length; i++) {
				String propValue = tokens[i].replace("\"", "");
				if(propValue.length() > 0) {
					mol.setProperty("ROOT:SALES(" + countMap.get(tokens[0]) + "):" + tags[i], propValue);
				}
			}
			countMap.put(tokens[0], Integer.valueOf(countMap.get(tokens[0]) + 1));
		}
		reader.close();
		resetCountMap(countMap);
		
		MolExporter molWriter = new MolExporter(new FileOutputStream(args[0]), "rdf");
		Set<String> molKeys = molMap.keySet();
		int[] molList = new int[molKeys.size()];
		int i = 0;
		for(String id : molKeys) {
			molList[i++] = Integer.parseInt(id);
		}
		java.util.Arrays.sort(molList);
		for(i = 0; i < molList.length; i++) {
			mol = molMap.get(Integer.toString(molList[i]));
			molWriter.write(mol);
		}
		molWriter.close();
	}
	
	public static void resetCountMap(HashMap<String, Integer> countMap) {
		for(Map.Entry<String, Integer> entry : countMap.entrySet()) {
			entry.setValue(Integer.valueOf(1));
		}
	}
	
}
