package com.sunset;

import java.io.IOException;
import java.util.HashMap;

import chemaxon.formats.MolExporter;
import chemaxon.struc.Molecule;

public class Fragment {
	
	private Molecule frag;
	private HashMap<String, Object> props = new HashMap<String, Object>();
	
	public Fragment(Molecule mol, int[] atomIndexes) {
		frag = mol.cloneMolecule();
		for(int n = frag.getAtomCount() - 1; n >= 0; n--) {
			if(Arrays.indexOf(n, atomIndexes) == -1) {
				frag.removeAtom(n);
			}
		}
		frag.pack();
	}
	
	public void addProp(String name, Object value) {
		props.put(name, value);
	}
	
	public Object getProp(String name) {
		return props.get(name);
	}
	
	public String toString(String fmt) {
		try {
			return MolExporter.exportToFormat(frag, fmt);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}
	
	public Molecule getMolecule() {
		return frag;
	}
	
}
