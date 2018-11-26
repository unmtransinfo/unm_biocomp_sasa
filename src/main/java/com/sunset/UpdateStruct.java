/**
 * 
 */
package com.sunset;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.io.MPropHandler;
import chemaxon.struc.Molecule;

/**
 * @author Oleg Ursu
 * update SDF file structure
 */
public class UpdateStruct {
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if(args.length != 4) {
			System.out.println("Usage: java UpdateStruct correct.sdf to_update.sdf tag_name output.sdf");
			return;
		}
		HashMap<String, Molecule> molMap = new HashMap<String, Molecule>();
		MolImporter reader = new MolImporter(args[0]);
		Molecule mol;
		while((mol = reader.read()) != null) {
			molMap.put(MPropHandler.convertToString(mol.properties(), args[2]), mol);
		}
		reader.close();
		reader = new MolImporter(args[1]);
		PrintWriter writer = new PrintWriter(args[3]);
		while((mol = reader.read()) != null) {
			if(molMap.containsKey(MPropHandler.convertToString(mol.properties(), args[2]))) {
				Molecule correctMol = molMap.get(mol.getPropertyObject(args[2]));
				for(String propKey : mol.properties().getKeys()) {
					correctMol.setProperty(propKey, MPropHandler.convertToString(mol.properties(), propKey));
				}
				writer.print(MolExporter.exportToFormat(correctMol, "sdf:-a-H"));
			} else {
				writer.print(MolExporter.exportToFormat(mol, "sdf:-a-H"));
			}
		}
		writer.close();
		reader.close();
	}

}
