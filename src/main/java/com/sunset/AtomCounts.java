package com.sunset;

import java.util.HashSet;

import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;

public class AtomCounts {
	
	private static final int[] counts = new int[105];
	private static final HashSet<Integer> nonMetal = new HashSet<Integer>();
	static {
		nonMetal.add(1);
		nonMetal.add(2);
		nonMetal.add(5);
		nonMetal.add(6);
		nonMetal.add(7);
		nonMetal.add(8);
		nonMetal.add(9);
		nonMetal.add(10);
		nonMetal.add(14);
		nonMetal.add(15);
		nonMetal.add(16);
		nonMetal.add(17);
		nonMetal.add(18);
		nonMetal.add(33);
		nonMetal.add(34);
		nonMetal.add(35);
		nonMetal.add(36);
		nonMetal.add(52);
		nonMetal.add(53);
		nonMetal.add(54);
		nonMetal.add(85);
		nonMetal.add(86);
	}
	
	private static final void reset() {
		for(int i = 0; i < counts.length; i++) {
			counts[i] = 0;
		}
	}

	public static final void elementCounts(Molecule mol) {
		reset();
		MolAtom[] atoms = mol.getAtomArray();
		for(MolAtom atom : atoms) {
			counts[atom.getAtno()]++;
		}
	}
	
	public static final int sumHetero() {
		return counts[7] + counts[8] + counts[15] + counts[16]; 
	}
	
	public static final int nonPolAtoms() {
		int pAtoms = sumHetero();
		int npAtoms = counts[6]+counts[9]+counts[17]+counts[35]+counts[53]; 
		if(npAtoms > pAtoms) {
			return npAtoms - pAtoms;
		} else {
			return 0;
		}
	}
	
	public static final int sumCarbons() {
		return counts[6];
	}

	public static final int sumMetal() {
		int metals = 0;
		for(int i = 1; i < counts.length; i++) {
			if(!nonMetal.contains(i) && counts[i] > 0) {
				metals++;
			}
		}
		return metals;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
	}

}
