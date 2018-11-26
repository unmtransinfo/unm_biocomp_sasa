/**
 * 
 */
package com.sunset;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.io.MPropHandler;
import chemaxon.struc.Molecule;
import chemaxon.struc.prop.MHashProp;
import chemaxon.struc.prop.MListProp;

/**
 * Split RDF input file in structure file SDF, ACT.LIST file and KW file
 * @author oleg
 *
 */
public class RDFSplit2 {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if(args.length != 1) {
			System.out.println("Usage: java RDFSplit RDF_file");
			return;
		}
		BufferedReader reader = new BufferedReader(new FileReader("wbpk_defs/drug_target.cols.txt"));
		ArrayList<String> drugTarget = new ArrayList<String>();
		String line;
		while((line = reader.readLine()) != null) {
			drugTarget.add(line);
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader("wbpk_defs/herg_bioact.cols.txt"));
		ArrayList<String> hergBioAct = new ArrayList<String>();
		while((line = reader.readLine()) != null) {
			hergBioAct.add(line);
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader("wbpk_defs/phase1_metab_enzyme.cols.txt"));
		ArrayList<String> phase1MetabEnzyme = new ArrayList<String>();
		while((line = reader.readLine()) != null) {
			phase1MetabEnzyme.add(line);
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader("wbpk_defs/wombat.cols.txt"));
		ArrayList<String> wombatActList = new ArrayList<String>();
		while((line = reader.readLine()) != null) {
			wombatActList.add(line);
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader("wbpk_defs/mixture.cols.txt"));
		ArrayList<String> mixture = new ArrayList<String>();
		while((line = reader.readLine()) != null) {
			mixture.add(line);
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader("wbpk_defs/mixture_ld50.cols.txt"));
		ArrayList<String> mixtureLD50 = new ArrayList<String>();
		while((line = reader.readLine()) != null) {
			mixtureLD50.add(line);
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader("wbpk_defs/mixture_water_sol.cols.txt"));
		ArrayList<String> mixtureWaterSol = new ArrayList<String>();
		while((line = reader.readLine()) != null) {
			mixtureWaterSol.add(line);
		}
		reader.close();
		
		PrintWriter drugTargetWriter = new PrintWriter("DRUG_TARGET.tab");
		drugTargetWriter.print("ID");
		for(int i = 0; i < drugTarget.size(); i++) {
			drugTargetWriter.print("\t" + drugTarget.get(i));
		}
		drugTargetWriter.println();
		
		PrintWriter hergBioActWriter = new PrintWriter("HERG_BIOACT.tab");
		hergBioActWriter.print("ID");
		for(int i = 0; i < hergBioAct.size(); i++) {
			hergBioActWriter.print("\t" + hergBioAct.get(i));
		}
		hergBioActWriter.println();
		
		
		
		PrintWriter phase1MetabEnzymeWriter = new PrintWriter("PHASE1_METAB_ENZYME.tab");
		phase1MetabEnzymeWriter.print("ID");
		for(int i = 0; i < phase1MetabEnzyme.size(); i++) {
			phase1MetabEnzymeWriter.print("\t" + phase1MetabEnzyme.get(i));
		}
		phase1MetabEnzymeWriter.println();
		
		
		PrintWriter wombatActListWriter = new PrintWriter("WOMBAT_ACTLIST.tab");
		wombatActListWriter.print("ID");
		for(int i = 0; i < wombatActList.size(); i++) {
			wombatActListWriter.print("\t" + wombatActList.get(i));
		}
		wombatActListWriter.println();
		
		
		PrintWriter mixtureWriter = new PrintWriter("MIXTURE.tab");
		mixtureWriter.print("ID\tWOMBAT_PK_MIXTURE_ID");
		for(int i = 0; i < mixture.size(); i++) {
			mixtureWriter.print("\t" + mixture.get(i));
		}
		mixtureWriter.println();
		
		
		PrintWriter mixtureLD50Writer = new PrintWriter("MIXTURE_LD50.tab");
		mixtureLD50Writer.print("WOMBAT_PK_MIXTURE_ID");
		for(int i = 0; i < mixtureLD50.size(); i++) {
			mixtureLD50Writer.print("\t" + mixtureLD50.get(i));
		}
		mixtureLD50Writer.println();
		
		
		
		PrintWriter mixtureWaterSolWriter = new PrintWriter("MIXTURE_WATER_SOL.tab");
		mixtureWaterSolWriter.print("WOMBAT_PK_MIXTURE_ID");
		for(int i = 0; i < mixtureWaterSol.size(); i++) {
			mixtureWaterSolWriter.print("\t" + mixtureWaterSol.get(i));
		}
		mixtureWaterSolWriter.println();
		
		
		Molecule mol;
		PrintWriter molWriter = new PrintWriter("WOMBAT-PK.sdf");
		int c = 0;
		int id = 0;
		int mixtureId = 0;
		MolImporter mReader = new MolImporter(new FileInputStream(args[0]), "rdf");
		while((mol = mReader.read()) != null) {
			c++;
			id = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "ROOT:ID"));
			Molecule mol2 = null;
			mol.properties().hierarchize();
			try {
				mol2 = mol.cloneMolecule();
				mol2.properties().clear();
			} catch (Exception e) {
				System.err.println("Can't read/parse record " + id);
				e.printStackTrace();
				continue;
			}
			if(mol.properties().get("ROOT") instanceof MHashProp) {
				MHashProp rootProps = (MHashProp)mol.properties().get("ROOT");
				for(int i = 0; i < rootProps.size(); i++) {
					if(rootProps.getKey(i).equals("DRUG_TARGET")) {
						MListProp subProps = (MListProp)rootProps.get(i);
						for(int k = 0; k < subProps.size(); k++) {
							drugTargetWriter.print(id);
							MHashProp prop = (MHashProp)subProps.get(k);
							for(int j = 0; j < drugTarget.size(); j++) {
								if(prop.get(drugTarget.get(j)) != null) {
									String str = prop.get(drugTarget.get(j)).getPropValue().toString();
									if(str.contains("\n")) {
										str = str.replace("\n", "");
									}
									drugTargetWriter.print("\t" + str);
								} else {
									drugTargetWriter.print("\t");
								}
							}
							drugTargetWriter.println();
						}
						continue;
					}
					
					if(rootProps.getKey(i).equals("HERG_BIOACT")) {
						MListProp subProps = (MListProp)rootProps.get(i);
						for(int k = 0; k < subProps.size(); k++) {
							
							hergBioActWriter.print(id);
							MHashProp prop = (MHashProp)subProps.get(k);
							for(int j = 0; j < hergBioAct.size(); j++) {
								if(prop.get(hergBioAct.get(j)) != null) {
									String str = prop.get(hergBioAct.get(j)).getPropValue().toString();
									if(str.contains("\n")) {
										str = str.replace("\n", "");
									}
									hergBioActWriter.print("\t" + str);
								} else {
									hergBioActWriter.print("\t");
								}
							}
							hergBioActWriter.println();
						}
						continue;
					}
					
					if(rootProps.getKey(i).equals("PHASE1_METAB_ENZYME")) {
						MListProp subProps = (MListProp)rootProps.get(i);
						for(int k = 0; k < subProps.size(); k++) {
							
							phase1MetabEnzymeWriter.print(id);
							MHashProp prop = (MHashProp)subProps.get(k);
							for(int j = 0; j < phase1MetabEnzyme.size(); j++) {
								if(prop.get(phase1MetabEnzyme.get(j)) != null) {
									String str = prop.get(phase1MetabEnzyme.get(j)).getPropValue().toString();
									if(str.contains("\n")) {
										str = str.replace("\n", "");
									}
									phase1MetabEnzymeWriter.print("\t" + str);
								} else {
									phase1MetabEnzymeWriter.print("\t");
								}
							}
							phase1MetabEnzymeWriter.println();
						}
						continue;
					}
					
					if(rootProps.getKey(i).equals("WOMBAT")) {
						MListProp subProps = (MListProp)rootProps.get(i);
						for(int k = 0; k < subProps.size(); k++) {
							
							wombatActListWriter.print(id);
							MHashProp prop = (MHashProp)subProps.get(k);
							for(int j = 0; j < wombatActList.size(); j++) {
								if(prop.get(wombatActList.get(j)) != null) {
									String str = prop.get(wombatActList.get(j)).getPropValue().toString();
									if(str.contains("\n")) {
										str = str.replace("\n", "");
									}
									wombatActListWriter.print("\t" + str);
								} else {
									wombatActListWriter.print("\t");
								}
							}
							wombatActListWriter.println();
						}
						continue;
					}
					
					if(rootProps.getKey(i).equals("MIXTURE")) {
						MListProp subProps = (MListProp)rootProps.get(i);
						for(int k = 0; k < subProps.size(); k++) {
							mixtureId++;
							MHashProp prop = (MHashProp)subProps.get(k);
							mixtureWriter.print(id + "\t" + mixtureId);
							for(int j = 0; j < mixture.size(); j++) {
								if(prop.get(mixture.get(j)) != null) {
									String str = prop.get(mixture.get(j)).getPropValue().toString();
									if(str.contains("\n")) {
										str = str.replace("\n", "");
									}
									mixtureWriter.print("\t" + str);
								} else {
									mixtureWriter.print("\t");
								}
							}
							mixtureWriter.println();
							if(prop.get("WATER_SOLUBILITY") != null) {
								MListProp subSubProp = (MListProp)prop.get("WATER_SOLUBILITY");
								for(int w = 0; w < subSubProp.size(); w++) {
									MHashProp prop2 = (MHashProp)subSubProp.get(w);
									mixtureWaterSolWriter.print(mixtureId);
									for(int u = 0; u < mixtureWaterSol.size(); u++) {
										if(prop2.get(mixtureWaterSol.get(u)) != null) {
											String str = prop2.get(mixtureWaterSol.get(u)).getPropValue().toString();
											if(str.contains("\n")) {
												str = str.replace("\n", "");
											}
											mixtureWaterSolWriter.print("\t" + str);
										} else {
											mixtureWaterSolWriter.print("\t");
										}
									}
									mixtureWaterSolWriter.println();
								}
							}
							if(prop.get("LD50") != null) {
								MListProp subSubProp = (MListProp)prop.get("LD50");
								for(int w = 0; w < subSubProp.size(); w++) {
									MHashProp prop2 = (MHashProp)subSubProp.get(w);
									mixtureLD50Writer.print(mixtureId);
									for(int u = 0; u < mixtureLD50.size(); u++) {
										if(prop2.get(mixtureLD50.get(u)) != null) {
											String str = prop2.get(mixtureLD50.get(u)).getPropValue().toString();
											if(str.contains("\n")) {
												str = str.replace("\n", "");
											}
											mixtureLD50Writer.print("\t" + str);
										} else {
											mixtureLD50Writer.print("\t");
										}
									}
									mixtureLD50Writer.println();
								}
							}
						}
						continue;
					}
					if(rootProps.getKey(i).equals("MOLNAME")) {
						mol2.setProperty(rootProps.getKey(i), rootProps.get(i).getPropValue().toString().toLowerCase());
					} else {
						mol2.setProperty(rootProps.getKey(i), rootProps.get(i).getPropValue().toString());
					}
				}
			}
			molWriter.print(MolExporter.exportToFormat(mol2, "sdf:-a-H"));
			if(c % 100 == 0) {
				System.out.println(c + " records processed");
			}
		}
		
		drugTargetWriter.close();
		hergBioActWriter.close();
		phase1MetabEnzymeWriter.close();
		wombatActListWriter.close();
		mixtureWriter.close();
		mixtureLD50Writer.close();
		mixtureWaterSolWriter.close();
		mReader.close();
		molWriter.close();		
	}

}
