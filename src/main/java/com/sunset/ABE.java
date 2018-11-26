/**
 * 
 */
package com.sunset;

import java.io.IOException;

import chemaxon.calculations.hydrogenize.Hydrogenize;
import chemaxon.formats.MolImporter;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.Molecule;

/**
 * Computes Andrew binding energy
 * @author Oleg Ursu
 *
 */
public class ABE {
	
	public static final Pattern[] patterns = new Pattern[] {
		new Pattern("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"),
		new Pattern("[$([OX2H1])!$(*a)!$(**=*)!$(*[#7]),$([O-1X1H0][C!H0])]"),
		new Pattern("[$([OX2H1]c),$([O-1X1H0]c)]"),
		new Pattern("[NH2D1!$(**=*)]-!@[#6!D1]"),
		new Pattern("[CX4v4]"),
		new Pattern("[#6X3v4+0,#6X3v3+1,#6X3v3-1]"),
		new Pattern("[#6X2v4+0,#6X1v3-1]"),
		new Pattern("[$([#6]=O)]"),
		new Pattern("[$([#8D1][#6]=[#8D1]),$([OH]C=O)]"),
		new Pattern("[$([#8X2H0][#6]=[#8D1])]"),
		new Pattern("[#7+1!$(*(=*)~*),$(C(N)=N)!$(C~a)!$(C=C)!$(*[#8]),$([N!R]cnc(N)n)]"),
		new Pattern("[#7X3!$(*a)!$(*-*=*)!$(*[#8])]"),
		new Pattern("[#7]"),
		new Pattern("[$([#15](~[#8])(~[#8])(~[#8])~[#8])]"),
		new Pattern("[F,Cl,Br,I]"),
		new Pattern("[#8,#16]"),
		new Pattern("[$([N!D1]-!@[C!$(*a)]=[N]),$([N!D1]-!@[C]=[O])]")
	};
	
	public static final double calcABE(Molecule mol) throws SearchException {
		double abe = -14.0d;
		int[] m = new int[patterns.length];
		int[][] matches = null;
		for(int i = 0; i < patterns.length; i++) {
			patterns[i].pattern.setTarget(mol);
			matches = patterns[i].pattern.findAll();
			if(matches != null) {
				m[i] = matches.length;
			}
		}
		// corrections
		// protonated nitrogen
		if(m[11] >= 1) {
			m[10]++;
		}
		// neutral nitrogen
		if(m[10] > 0) {
			m[12] -= m[10];
		}
		abe += -0.7 * (m[0] + m[1] + m[2] + m[3] - m[16]) + 0.7 * (m[5] - m[7] + m[6]) + 0.8 * m[4] + 11.5 * m[10] +
			1.2 * m[12] + 8.2 * m[8] + 10.0 * m[13] + 2.5 * (m[1] + m[2]) + 3.4 * (m[7] + m[9]) + 1.1 * (m[15] + m[9]) +
			1.3 * m[14];
		return abe;
	} 

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException, SearchException {
		Molecule mol = MolImporter.importMol("CCCC[C@H]1CN(C(=O)C(=O)N1Cc1cncn1Cc1ccc(cc1)C#N)c1cccc(C)c1C", "smiles");
		Hydrogenize.removeHAtoms(mol);
		System.out.println(ABE.calcABE(mol));
	}

}
