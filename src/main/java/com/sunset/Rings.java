package com.sunset;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import chemaxon.formats.MolImporter;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;


public class Rings {

	public static final Pattern[] RINGS = new Pattern[] {
		new Pattern("*~1~*~*~1"),
		new Pattern("*~1~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~*~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~*~*~*~*~*~1")
	};
	
	private static int[][] sssr = null;
	private static int[][] sssrBonds = null;
	private static Molecule mol;
	
	public static void setMolecule(Molecule mol) {
		Rings.mol = mol;
		sssr = mol.getSSSR();
		sssrBonds = new int[sssr.length][];
		for(int r = 0; r < sssr.length; r++) {
			sssrBonds[r] = getRingBonds(sssr[r]);
		}
	}
	
	// private static TopologyAnalyser tp = new TopologyAnalyser();
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException, SearchException {
		Molecule mol = MolImporter.importMol("O=C(CC1CC(=O)NC(CC2=CNC3=C2C=CC=C3)C(=O)NC(CC2=CC=CC=C2)C(=O)NC(CC2=CC=CC=C2)CNC1=O)N1CCOCC1", "smiles");
		mol.aromatize(MoleculeGraph.AROM_GENERAL);
		setMolecule(mol);
		System.out.println(heteroAliphaticRings() != null ? heteroAliphaticRings().length : 0);
		System.out.println(heteroAromaticRings() != null ? heteroAromaticRings().length : 0);
	}
	
	public static int[][] getSSSR() {
		return sssr;
	}
	
	public static int[] getRingBonds(int[] ringAtoms) {
		int[][] btab = mol.getBondTable().getMatrixArray();
		HashSet<Integer> ringBonds = new HashSet<Integer>();
		for(int i = 0; i < ringAtoms.length; i++) {
			for(int j = i + 1; j < ringAtoms.length; j++) {
				if(btab[ringAtoms[i]][ringAtoms[j]] == -1) {
					continue;
				} else {
					ringBonds.add(btab[ringAtoms[i]][ringAtoms[j]]);
				}
			}
		}
		int[] ans = new int[ringBonds.size()];
		int c = 0;
		for(Integer ringBond : ringBonds) {
			ans[c++] = ringBond;
		}
		return ans;
	}
	
	
	public static int[][] heteroAliphaticRings() {
		// tp.setMolecule(mol, MoleculeGraph.AROM_GENERAL);
		ArrayList<Integer> heteroAliphaticRingsIdx = new ArrayList<Integer>();
		int[][] aliphaticRings = aliphaticRings();
		if(aliphaticRings != null) {
			for(int j = 0; j < aliphaticRings.length; j++) {
				for(int i : aliphaticRings[j]) {
					if(mol.getAtom(i).getAtno() == 7 || mol.getAtom(i).getAtno() == 8 ||
							mol.getAtom(i).getAtno() == 15 || mol.getAtom(i).getAtno() == 16) {
							heteroAliphaticRingsIdx.add(j);
							break;
					}
				}

			}			
		} else {
			return null;
		}
		int[][] heteroAliphaticRings = new int[heteroAliphaticRingsIdx.size()][];
		for(int i = 0; i < heteroAliphaticRings.length; i++) {
			heteroAliphaticRings[i] = aliphaticRings[heteroAliphaticRingsIdx.get(i)];
		}
		return heteroAliphaticRings;
	}
	
	public static int[][] heteroAromaticRings() {
		// tp.setMolecule(mol, MoleculeGraph.AROM_GENERAL);
		ArrayList<Integer> heteroAromaticRingsIdx = new ArrayList<Integer>();
		int[][] aromaticRings = aromaticRings();
		if(aromaticRings != null) {
			for(int j = 0; j < aromaticRings.length; j++) {
				for(int i : aromaticRings[j]) {
					if(mol.getAtom(i).getAtno() == 7 || mol.getAtom(i).getAtno() == 8 ||
							mol.getAtom(i).getAtno() == 15 || mol.getAtom(i).getAtno() == 16) {
							heteroAromaticRingsIdx.add(j);
							break;
					}
				}

			}			
		} else {
			return null;
		}
		int[][] heteroAromaticRings = new int[heteroAromaticRingsIdx.size()][];
		for(int i = 0; i < heteroAromaticRings.length; i++) {
			heteroAromaticRings[i] = aromaticRings[heteroAromaticRingsIdx.get(i)];
		}
		return heteroAromaticRings;
	}
	
	public static int[][] aliphaticRings() {
		ArrayList<int[]> aliphaticRings = new ArrayList<int[]>();
		for(int r = 0; r < sssrBonds.length; r++) {
			int[] ringBonds = sssrBonds[r];
			if(!isAromaticRing(ringBonds)) {
				aliphaticRings.add(sssr[r]);
			}
		}
		if(aliphaticRings.size() == 0) {
			return null;
		}
		int[][] ans = new int[aliphaticRings.size()][];
		for(int i = 0; i < aliphaticRings.size(); i++) {
			ans[i] = aliphaticRings.get(i);
		}
		return ans;
	}
	
	public static int[][] aromaticRings() {
		ArrayList<int[]> aromaticRings = new ArrayList<int[]>();
		for(int r = 0; r < sssrBonds.length; r++) {
			int[] ringBonds = sssrBonds[r];
			if(isAromaticRing(ringBonds)) {
				aromaticRings.add(sssr[r]);
			}
		}
		if(aromaticRings.size() == 0) {
			return null;
		}
		int[][] ans = new int[aromaticRings.size()][];
		for(int i = 0; i < aromaticRings.size(); i++) {
			ans[i] = aromaticRings.get(i);
		}
		return ans;
	}
	
	public static boolean isAromaticRing(int[] ringBonds) {
		for(int b = 0; b < ringBonds.length; b++) {
			MolBond bond = mol.getBond(ringBonds[b]);
			if(bond.getType() != 4) {
				return false;
			} 
		}
		return true;
	}
	
}
