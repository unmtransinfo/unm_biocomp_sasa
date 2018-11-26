/**
 * 
 */
package com.sunset;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import chemaxon.sss.search.SearchException;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;

/**
 * Computes number of rotatable and rigid bonds in a molecule
 * @author Oleg Ursu
 *
 */
public class RotRigBonds {
	
	public static final Pattern AMIDE = new Pattern("[O,N][C,N,S,P]=[O,N]");
	
	
	public static final void rotBonds(Molecule mol, int[] rt) throws SearchException {
		rt[0] = rt[1] = 0;
		HashSet<Integer> termBonds = findTerminalBonds(mol);
		int[][] btab = mol.getBondTable().getMatrixArray();
		HashSet<Integer> rigidBonds = new HashSet<Integer>();
		int[][] ringBonds = getRingBonds(mol);
		int bondCount = 0;
		int size = 0;
		int contrib = 0;
		
		/* generate a list of rigid amide bonds              */
		/* (other rigid bonds can be included if necessary)  */
		AMIDE.pattern.setTarget(mol);
		int[][] matches = null;
		matches = AMIDE.pattern.findAll();
		if(matches != null) {
			for(int[] match : matches) {
				rigidBonds.add(btab[match[0]][match[1]]);
				rigidBonds.add(btab[match[1]][match[2]]);
			}
		}
		
		/* consider all acyclic bonds first */
		for(int i = 0; i < mol.getBondCount(); i++) {
			/* eliminate all non-single bonds... */
			MolBond bond = mol.getBond(i);
			if(bond.getType() != 1) {
				continue;
			}
			/* ... terminal bonds ... */
			if(termBonds.contains(i)) {
				continue;
			}
			/* ... and bonds which are incorporated in a ring system ... */
			if(mol.isRingBond(i)) {
				continue;
			}
			/* ... and not an amide bond */
			if(rigidBonds.contains(i)) {
				continue;
			}
			rt[0]++;
		}
		
		/* consider all rings */
		for(int[] ring : ringBonds) {
			bondCount = 0;
			size = ring.length;
			for(int bond : ring) {
				if(mol.getBond(bond).getType() != 1) {
					bondCount++;
					continue;
				}
				if(countBondRings(ringBonds, bond) > 1) {
					bondCount++;
					continue;
				}
				if(rigidBonds.contains(bond)) {
					bondCount++;
				}
			}
			contrib = size - 4 - bondCount;
			if(contrib < 0) {
				contrib = 0;
			}
			rt[0] += contrib;
		}
		rt[1] = termBonds.size();
	}
	
	public static final HashSet<Integer> findTerminalBonds(Molecule mol) {
		HashSet<Integer> tBonds = new HashSet<Integer>();
		for(int i = 0; i < mol.getBondCount(); i++) {
			if(mol.getBond(i).getType() != 1) {
				continue;
			}
			if(mol.getBond(i).getAtom1().getBondCount() == 1 ||
			   mol.getBond(i).getAtom2().getBondCount() == 1 ||
			   mol.getBond(i).getAtom1().getBondCount() - mol.getBond(i).getAtom1().getExplicitHcount() < 2 ||
			   mol.getBond(i).getAtom2().getBondCount() - mol.getBond(i).getAtom2().getExplicitHcount() < 2) {
				tBonds.add(i);
			}
		}
		return tBonds;
	}
	
	public static final int[][] getRingBonds(Molecule mol) {
		int[][] aRings = mol.getSSSR();
		int[][] bRings = new int[aRings.length][];
		int[][] bonds = mol.getBondTable().getMatrixArray();
		
		for(int r = 0; r < aRings.length; r++) {
			ArrayList<Integer> bRingList = new ArrayList<Integer>();
			int i1 = 0;
			while(bRingList.size() != aRings[r].length) {
				for(int k = 0; k < aRings[r].length; k++) {
					if(k == i1) {
						continue;
					}
					int bond = bonds[aRings[r][i1]][aRings[r][k]]; 
					if(bond != -1 && bRingList.indexOf(bond) == -1) {
						bRingList.add(bond);
						i1 = k;
						break;
					}
				}
			}
			int[] bRing = new int[bRingList.size()];
			for(int i = 0; i < bRing.length; i++) {
				bRing[i] = bRingList.get(i);
			}
			bRings[r] = bRing;
		}
		return bRings;
	}
	
	public static final int countBondRings(int[][] ringBonds, int bond) {
		int c = 0;
		for(int[] ring : ringBonds) {
			for(int i : ring) {
				if(i == bond) {
					c++;
				}
			}
		}
		return c;
	}
	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException, SearchException {
		MolHandler mh = new MolHandler("O=C1NC2=C(C=CC=C2N(=O)=O)C2=C1C=CC=C2");
		int[] rt = new int[2];
		rotBonds(mh.getMolecule(), rt);
		System.out.println(Arrays.toString(rt));
		System.out.println(mh.getMolecule().getBondCount() - rt[0] - rt[1]);
	}

}
