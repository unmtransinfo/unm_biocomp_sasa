package com.sunset;

import chemaxon.calculations.hydrogenize.Hydrogenize;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;

/**
 * @author Oleg Ursu
 * Compute Atomic and Bond Contributions of van der Waals volume (VABC)
 * according to method described in http://dx.doi.org/10.1021/jo034808o
 */
public class VABC {
	
	/**
	 * calculate van der Waals volume based on atom and bond contributions
	 * @param mol input structure
	 * @return volume (A) or 0 if volume can not be computed
	 */
	
	public static double calculate(Molecule mol) {
		double volume = 0.0d;
		for(MolAtom atom : mol.getAtomArray()) {
			if(VDW.ofNumber(atom.getAtno()) == null) {
				return 0.0d;
			}
			volume += VDW.ofNumber(atom.getAtno()).volume();
		}
		Rings.setMolecule(mol);
		int[][] aromaticRings = Rings.aromaticRings();
		int[][] aliphaticRings = Rings.aliphaticRings();
		volume = volume - 5.92d*mol.getBondCount() - 14.7d*(aromaticRings == null ? 0 : aromaticRings.length) - 3.8d*(aliphaticRings == null ? 0 : aliphaticRings.length);
		return volume;
	}

	/**
	 * @param args
	 * @throws MolFormatException 
	 */
	public static void main(String[] args) throws MolFormatException {
		// TODO Auto-generated method stub
		Molecule mol = MolImporter.importMol("Nc1nccs1", "smiles");
		mol.aromatize(MoleculeGraph.AROM_GENERAL);
		Hydrogenize.addHAtoms(mol);
		System.out.println(mol.getBondCount());
		System.out.println(VABC.calculate(mol));
	}

}
