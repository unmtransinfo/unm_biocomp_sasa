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
import java.util.Arrays;

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
public class RDFSplit {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if(args.length != 1) {
			System.out.println("Usage: java RDFSplit RDF_file");
			return;
		}
		/*boolean select = args[2].equals("null") ? false : true;
		String selectProp = "";
		HashSet<String> selectProps = new HashSet<String>();
		if(select) {
			BufferedReader selPropReader = new BufferedReader(new FileReader(args[2]));
			String line;
			selectProp = selPropReader.readLine();
			while((line = selPropReader.readLine()) != null) {
				selectProps.add(line);
			}
			selPropReader.close();
		}*/
		ArrayList<String> actListFields = new ArrayList<String>();
		BufferedReader reader = new BufferedReader(new FileReader("wb_defs/act.list.fields.txt"));
		String line;
		while((line = reader.readLine()) != null) {
			actListFields.add(line);
		}
		reader.close();

		reader = new BufferedReader(new FileReader("wb_defs/vivo.act.list.fields.txt"));
		ArrayList<String> vivoActListFields = new ArrayList<String>();
		while((line = reader.readLine()) != null) {
			vivoActListFields.add(line);
		}
		reader.close();

		reader = new BufferedReader(new FileReader("wb_defs/admet.act.list.fields.txt"));
		ArrayList<String> admetActListFields = new ArrayList<String>();
		while((line = reader.readLine()) != null) {
			admetActListFields.add(line);
		}
		reader.close();

		PrintWriter actWriter = new PrintWriter("ACT.LIST.TAB");
		actWriter.print("SMDL.ID");
		for(int i = 0; i < actListFields.size(); i++) {
			actWriter.print("\t" + actListFields.get(i));
		}
		actWriter.println();

		PrintWriter vivoActWriter = new PrintWriter("VIVO.ACT.LIST.TAB");
		vivoActWriter.print("SMDL.ID");
		for(int i = 0; i < vivoActListFields.size(); i++) {
			vivoActWriter.print("\t" + vivoActListFields.get(i));
		}
		vivoActWriter.println();

		PrintWriter admetActWriter = new PrintWriter("ADMET.ACT.LIST.TAB");
		admetActWriter.print("SMDL.ID");
		for(int i = 0; i < admetActListFields.size(); i++) {
			admetActWriter.print("\t" + admetActListFields.get(i));
		}
		admetActWriter.println();

		MolImporter mReader = new MolImporter(new FileInputStream(args[0]), "rdf");
		Molecule mol;
		PrintWriter molWriter = new PrintWriter("WOMBAT.sdf");
		PrintWriter kwWriter = new PrintWriter("KW.TAB");
		kwWriter.println("SMDL.ID\tKW");
		int c = 0;
		int smdlId = 0;
		String str = "";
		while((mol = mReader.read()) != null) {
			c++;
			/*if(select) {
				String value = mol.getProperty(selectProp);
				if(!selectProps.contains(value)) {
					continue;
				}
			}*/
			smdlId = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "ROOT:SMDL.ID"));
			Molecule mol2 = null;
			try {
				mol2 = mol.cloneMolecule();
				mol2.properties().clear();
			} catch (Exception e) {
				System.err.println("error parsing SMDL.ID record " + smdlId);
				System.err.println(Arrays.toString(e.getStackTrace()));
				System.exit(1);
			}
			mol.properties().hierarchize();
			if(mol.properties().get("ROOT") instanceof MHashProp) {
				MHashProp rootProps = (MHashProp)mol.properties().get("ROOT");
				for(int i = 0; i < rootProps.size(); i++) {
					if(rootProps.getKey(i).equals("ACT.LIST")) {
						MListProp actList = (MListProp)rootProps.get(i);
						for(int k = 0; k < actList.size(); k++) {
							if(actList.get(k) == null) {
								continue;
							}
							actWriter.print(rootProps.get("SMDL.ID").getPropValue());
							MHashProp act = (MHashProp)actList.get(k);
							for(int j = 0; j < actListFields.size(); j++) {
								if(act.get(actListFields.get(j)) != null) {
									str = act.get(actListFields.get(j)).getPropValue().toString();
									if(str.contains("\n")) {
										str = str.replace("\n", "");
									}
									str = str.replaceAll("\\s+", " ");
									actWriter.print("\t" + str);
								} else {
									actWriter.print("\t");
								}
							}
							actWriter.println();
						}
						continue;
					}
					if(rootProps.getKey(i).equals("ACT.LIST1")) {
						MListProp actList = (MListProp)rootProps.get(i);
						for(int k = 0; k < actList.size(); k++) {
							vivoActWriter.print(rootProps.get("SMDL.ID").getPropValue());
							MHashProp act = (MHashProp)actList.get(k);
							for(int j = 0; j < vivoActListFields.size(); j++) {
								if(act.get(vivoActListFields.get(j)) != null) {
									str = act.get(vivoActListFields.get(j)).getPropValue().toString();
									if(str.contains("\n")) {
										str = str.replace("\n", "");
									}
									str = str.replaceAll("\\s+", " ");
									vivoActWriter.print("\t" + str);
								} else {
									vivoActWriter.print("\t");
								}
							}
							vivoActWriter.println();
						}
						continue;
					}
					if(rootProps.getKey(i).equals("ACT.LIST2")) {
						MListProp actList = (MListProp)rootProps.get(i);
						for(int k = 0; k < actList.size(); k++) {
							admetActWriter.print(rootProps.get("SMDL.ID").getPropValue());
							MHashProp act = (MHashProp)actList.get(k);
							for(int j = 0; j < admetActListFields.size(); j++) {
								if(act.get(admetActListFields.get(j)) != null) {
									str = act.get(admetActListFields.get(j)).getPropValue().toString();
									if(str.contains("\n")) {
										str = str.replace("\n", "");
									}
									str = str.replaceAll("\\s+", " ");
									admetActWriter.print("\t" + str);
								} else {
									admetActWriter.print("\t");
								}
							}
							admetActWriter.println();
						}
						continue;
					}
					if(rootProps.getKey(i).equals("MOL.KW")) {
						MListProp kwList = (MListProp)rootProps.get(i);
						for(int k = 0; k < kwList.size(); k++) {
							kwWriter.print(rootProps.get("SMDL.ID").getPropValue());
							MHashProp kw = (MHashProp)kwList.get(k);
							kwWriter.println("\t" + kw.get("KW").getPropValue());
						}
						continue;
					}
					mol2.setProperty(rootProps.getKey(i), rootProps.get(i).getPropValue().toString());
				}
			}
			molWriter.print(MolExporter.exportToFormat(mol2, "sdf:-a-H"));
			if(c % 1000 == 0) {
				System.out.println(c + " records processed");
			}
		}
		mReader.close();
		actWriter.close();
		vivoActWriter.close();
		admetActWriter.close();
		molWriter.close();
		kwWriter.close();		
	}

}
