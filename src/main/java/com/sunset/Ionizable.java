/**
 * 
 */
package com.sunset;

import java.io.IOException;
import java.util.Arrays;

import chemaxon.sss.search.MolSearch;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;
import chemaxon.util.MolHandler;

/**
 * Computes pozitive/negative ionizable atoms in molecule adapted from SaSA 
 * @author Oleg Ursu
 */
public class Ionizable {

	private static final MolHandler mh = new MolHandler();

	private static final MolSearch ms = new MolSearch();
	
	static {
		mh.setQueryMode(true);
	}
	
	public static final void ionizableGroups(Molecule mol, int[] pn) throws IOException, SearchException {
		pn[0] = pn[1] = 0;
		int[][] matches = null;
		ms.setTarget(mol);
		
		/* negatively ionized acids */
		mh.setMolecule("[OH1;$(O[C,S,P](=O))]");
		ms.setQuery(mh.getMolecule());
		matches = ms.findAll();
		if(matches != null) {
			pn[1] += matches.length;
		}
		
		/* negatively charged atoms, ex those at semipolar bonds */
		mh.setMolecule("[-;!$([-]~[+])]");
		ms.setQuery(mh.getMolecule());
		matches = ms.findAll();
		if(matches != null) {
			pn[1] += matches.length;
		}
		
		/* positively ionized nitrogens */
		
		/* amines, ex. amides, nitro/nitroso, amidines, guanidines */
		mh.setMolecule("[N;+0;!$(NC=O);!$(NS=O);!$(N=CN);!$(NC=N);!$(NC=S);!$(N~O)]");
		ms.setQuery(mh.getMolecule());
		matches = ms.findAll();
		if(matches != null) {
			pn[0] += matches.length;
		}
		
		/* pyridine N, ex pyrrole, pyridinium salts */
		mh.setMolecule("[n;-0]1ccccc1");
		ms.setQuery(mh.getMolecule());
		matches = ms.findAll();
		if(matches != null) {
			pn[0] += matches.length;
		}
		
		/* guanidines */
		mh.setMolecule("[$([N;+0]C(=[N;+0])[N;+0])]");
		ms.setQuery(mh.getMolecule());
		matches = ms.findAll();
		if(matches != null) {
			pn[0] += matches.length / 2;
		}
		
		/* amidines, ex. guanidines */
		mh.setMolecule("[$([N;+0]C=[N;+0]);!$(NC(=N)N)]");
		ms.setQuery(mh.getMolecule());
		matches = ms.findAll();
		if(matches != null) {
			pn[0] += matches.length;
		}
		
		/* positively charged atoms, ex those at semipolar bonds */
		mh.setMolecule("[+;!$([+]~[-])]");
		ms.setQuery(mh.getMolecule());
		matches = ms.findAll();
		if(matches != null) {
			pn[0] += matches.length;
		}
	}	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException, SearchException {
		MolHandler mh1 = new MolHandler("CC(C)C[C@@H]1NC(=O)[C@@H](NC(=O)[C@H](CC2=CNC3=C2C=CC=C3)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CCCCN)NC(=O)[C@H](CCCNC(N)=N)NC1=O)C(C)(C)C");
		int[] pn = new int[2];
		mh1.aromatize(MoleculeGraph.AROM_GENERAL);
		ionizableGroups(mh1.getMolecule(), pn);
		System.out.println("Positive, Negative");
		System.out.println(Arrays.toString(pn));
		AtomCounts.elementCounts(mh1.getMolecule());
		System.out.println("sum hetero: ");
		System.out.println(AtomCounts.sumHetero());
		System.out.println("non_polar_atoms: ");
		System.out.println(AtomCounts.nonPolAtoms());
		System.out.println("ABE: ");
		System.out.println(ABE.calcABE(mh1.getMolecule()));
	}

}
