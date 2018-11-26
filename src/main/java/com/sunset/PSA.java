/**
 * 
 */
package com.sunset;

import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

import chemaxon.calculations.hydrogenize.Hydrogenize;
import chemaxon.formats.MolImporter;
import chemaxon.sss.search.MolSearch;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;
import chemaxon.util.MolHandler;


/**
 * @author Oleg Ursu
 *
 */
public class PSA {
	
	public static final double FLT_PREC = 1.0e-7d;

	public static final VDWType[] VDW = new VDWType[] {
			new VDWType("[#1]", Bool.UNKNOWN, 1.20d),
			new VDWType("[#6]", Bool.FALSE, 1.70d),
			new VDWType("[#7]", Bool.TRUE, 1.55d),
			new VDWType("[#8]", Bool.TRUE, 1.52d),
			new VDWType("[S,s]", Bool.FALSE, 1.80d),
			new VDWType("[P]", Bool.FALSE, 1.80d),
			new VDWType("[F]", Bool.FALSE, 1.40d),
			new VDWType("[Cl]", Bool.FALSE, 1.75d),
			new VDWType("[Br]", Bool.FALSE, 1.95d),
			new VDWType("[I]", Bool.FALSE, 2.10d)
		};
	
	public static final D12Type[] D12 = new D12Type[] {
		new D12Type("[#1]-[#6]", 1.11d),
		new D12Type("[#1]-[#7]", 1.01d),
		new D12Type("[#1]-O", 0.96d),
		new D12Type("N#C", 1.16d),
		new D12Type("C#C", 1.20d),
		new D12Type("N=N", 1.25d),
		new D12Type("[O,N]=[#6]", 1.25d),
		new D12Type("O=[#7]", 1.30d),
		new D12Type("S-[#1]", 1.30d),
		new D12Type("C=C", 1.35d),
		new D12Type("O-C", 1.35d),
		new D12Type("F-c", 1.35d),
		new D12Type("O-[#7]", 1.35d),
		new D12Type("n:n", 1.35d),
		new D12Type("[N,O,n]-c", 1.40d),
		new D12Type("F-C", 1.40d),
		new D12Type("c:[c,n,o]", 1.40d),
		new D12Type("n-N", 1.40d),
		new D12Type("o:n", 1.40d),
		new D12Type("N-[C,N]", 1.45d),
		new D12Type("S=O", 1.45d),
		new D12Type("c-c", 1.45d),
		new D12Type("n-C", 1.45d),
		new D12Type("C-[#6]", 1.54d),
		new D12Type("O-O", 1.50d),
		new D12Type("P-F", 1.50d),
		new D12Type("P=O", 1.50d),
		new D12Type("P-O", 1.55d),
		new D12Type("s:n", 1.55d),
		new D12Type("S-[O,F]", 1.60d),
		new D12Type("Cl-c", 1.70d),
		new D12Type("P-N", 1.70d),
		new D12Type("S-n", 1.70d),
		new D12Type("S=[#6]", 1.70d),
		new D12Type("s:c", 1.70d),
		new D12Type("Cl-[C,#7]", 1.75d),
		new D12Type("S-N", 1.75d),
		new D12Type("P-[#6]", 1.80d),
		new D12Type("S-[#6]", 1.80d),
		new D12Type("[#14]-C", 1.85d),
		new D12Type("Br-[C,c]", 1.90d),
		new D12Type("S=P", 1.95d),
		new D12Type("S-S", 2.00d),
		new D12Type("I-c", 2.05d),
		new D12Type("I-C", 2.15d)
	};
	
	public static final int getPSA(Molecule mol, double[] psa) throws SearchException, IOException {
		Hydrogenize.addHAtoms(mol);
		double[][] dm = new double[mol.getAtomCount()][mol.getAtomCount()];
		getDist_12(mol, dm);
		calcPSA(mol, dm, psa);
		return 0;
	}
	
	public static final int getPSA(Molecule mol, double[] psa, int[] atomIndexes) throws SearchException, IOException {
		Hydrogenize.addHAtoms(mol);
		double[][] dm = new double[mol.getAtomCount()][mol.getAtomCount()];
		getDist_12(mol, dm);
		calcPSA(mol, dm, psa, atomIndexes);
		return 0;
	}
	
	public static final int calcPSA(Molecule mol, double[][] dm, double[] psa, int[] atomIndexes) throws IOException, SearchException {
		psa[0] = psa[1] = 0.0d;
		double[] radius = new double[mol.getAtomCount()];
		double[] surface = new double[mol.getAtomCount()];
		Bool[] polar = new Bool[mol.getAtomCount()];
		double[] s12 = new double[2];
		for(int i = 0; i < VDW.length; i++) {
			VDW[i].pattern.setTarget(mol);
			int[][] matches = VDW[i].pattern.findAll();
			if(matches == null) {
				continue;
			}
			for(int[] match : matches) {
				int natom1 = match[0];
				radius[natom1] = VDW[i].radius;
				polar[natom1] = VDW[i].polar;
			}
		}
		// assign polarity type to Hs
		MolHandler mh = new MolHandler();
		mh.setQueryMode(true);
		mh.setMolecule("[#7,#8][#1]");
		MolSearch pattern = new MolSearch();
		pattern.setQuery(mh.getMolecule());
		pattern.setTarget(mol);
		int[][] matches = pattern.findAll();
		if(matches != null) {
			for(int[] match : matches) {
				int natom2 = match[1];
				polar[natom2] = Bool.TRUE;
			}
		}
		mh.setMolecule("[!#7&!#8][#1]");
		pattern.setQuery(mh.getMolecule());
		pattern.setTarget(mol);
		matches = pattern.findAll();
		if(matches != null) {
			for(int[] match : matches) {
				int natom2 = match[1];
				polar[natom2] = Bool.FALSE;
			}
		}
		// summing vdW radii for polar and non-polar atoms
		for(int i = 0; i < surface.length; i++) {
			surface[i] = 4.0d * Math.PI * (radius[i] * radius[i]);
		}
		
		// correction for bonded atoms overlap
		mh.setMolecule("[*,#1]~[*,#1]");
		pattern.setQuery(mh.getMolecule());
		pattern.setTarget(mol);
		matches = pattern.findAll();
		if(matches != null) {
			for(int[] match : matches) {
				int natom1 = match[0];
				int natom2 = match[1];
				if(Math.abs(radius[natom1]) > FLT_PREC && Math.abs(radius[natom2]) > FLT_PREC &&
						Math.abs(dm[natom1][natom2]) > FLT_PREC) {
					hiddenSurface12(radius[natom1],radius[natom2],dm[natom1][natom2],s12);
				    surface[natom1] -= s12[0];
				    surface[natom2] -= s12[1];
				}
			}
		}
		int[][] bondTable = mol.getBondTable().getMatrixArray();
		for(int i = 0; i < mol.getAtomCount(); i++) {
			if(surface[i] < 0.0D) {
				surface[i] = 0.0d;
			}
			if(com.sunset.Arrays.indexOf(i, atomIndexes) == -1 && mol.getAtom(i).getAtno() != 1) {
				continue;
			}
			if(mol.getAtom(i).getAtno() == 1) {
				boolean connected = false;
				for(int k : atomIndexes) {
					if(bondTable[i][k] != -1) {
						connected = true;
					}
				}
				if(!connected) {
					continue;
				}
			}
			// sum up polar/nonpolar atomic contributions
			if(polar[i] == Bool.TRUE) {
				psa[0] += surface[i];
			} else {
				psa[1] += surface[i];
			}
		}
		return 0;
	}
	
	public static final int calcPSA(Molecule mol, double[][] dm, double[] psa) throws IOException, SearchException {
		psa[0] = psa[1] = 0.0d;
		double[] radius = new double[mol.getAtomCount()];
		double[] surface = new double[mol.getAtomCount()];
		Bool[] polar = new Bool[mol.getAtomCount()];
		double[] s12 = new double[2];
		for(int i = 0; i < VDW.length; i++) {
			VDW[i].pattern.setTarget(mol);
			int[][] matches = VDW[i].pattern.findAll();
			if(matches == null) {
				continue;
			}
			for(int[] match : matches) {
				int natom1 = match[0];
				radius[natom1] = VDW[i].radius;
				polar[natom1] = VDW[i].polar;
			}
		}
		// assign polarity type to Hs
		MolHandler mh = new MolHandler();
		mh.setQueryMode(true);
		mh.setMolecule("[#7,#8][#1]");
		MolSearch pattern = new MolSearch();
		pattern.setQuery(mh.getMolecule());
		pattern.setTarget(mol);
		int[][] matches = pattern.findAll();
		if(matches != null) {
			for(int[] match : matches) {
				int natom2 = match[1];
				polar[natom2] = Bool.TRUE;
			}
		}
		mh.setMolecule("[!#7&!#8][#1]");
		pattern.setQuery(mh.getMolecule());
		pattern.setTarget(mol);
		matches = pattern.findAll();
		if(matches != null) {
			for(int[] match : matches) {
				int natom2 = match[1];
				polar[natom2] = Bool.FALSE;
			}
		}
		// summing vdW radii for polar and non-polar atoms
		for(int i = 0; i < surface.length; i++) {
			surface[i] = 4.0d * Math.PI * (radius[i] * radius[i]);
		}
		
		// correction for bonded atoms overlap
		mh.setMolecule("[*,#1]~[*,#1]");
		pattern.setQuery(mh.getMolecule());
		pattern.setTarget(mol);
		matches = pattern.findAll();
		if(matches != null) {
			for(int[] match : matches) {
				int natom1 = match[0];
				int natom2 = match[1];
				if(Math.abs(radius[natom1]) > FLT_PREC && Math.abs(radius[natom2]) > FLT_PREC &&
						Math.abs(dm[natom1][natom2]) > FLT_PREC) {
					hiddenSurface12(radius[natom1],radius[natom2],dm[natom1][natom2],s12);
				    surface[natom1] -= s12[0];
				    surface[natom2] -= s12[1];
				}
			}
		}
		for(int i = 0; i < mol.getAtomCount(); i++) {
			if(surface[i] < 0.0D) {
				surface[i] = 0.0d;
			}
			// sum up polar/nonpolar atomic contributions
			if(polar[i] == Bool.TRUE) {
				psa[0] += surface[i];
			} else {
				psa[1] += surface[i];
			}
		}
		return 0;
	}
	
	public static final void getDist_12(Molecule mol, double[][] dm) throws SearchException {
		for(int i = 0; i < D12.length; i++) {
			D12[i].pattern.setTarget(mol);
			int[][] matches = D12[i].pattern.findAll();
			if(matches == null) {
				continue;
			}
			for(int[] match : matches) {
				int natom1 = match[0];
				int natom2 = match[1];
				if(Math.abs(dm[natom1][natom2]) < FLT_PREC) {
					dm[natom1][natom2] = D12[i].dist;
					dm[natom2][natom1] = D12[i].dist;
				}
			}
		}
	}
	
	public static final int hiddenSurface12(double r1, double r2, double d, double[] s12) {
		s12[0] = s12[1] = 0.0d;
		if(r1 <= 0.0D || r2 <= 0.0D || d <= 0.0D) {
			return -1;
		}
		// spheres do not intersect
		if(d > r1 + r2) {
			return 0;
		}
		// one sphere is completely inside another one
		if(d < Math.abs(r1 - r2)) {
			if(r1 < r2) {
				s12[0] = 4. * Math.PI * r1 * r1;
				s12[1] = 0.0;
			} else {
				s12[0] = 0.0;
			    s12[1] = 4. * Math.PI * r2 * r2;
			}
			return 0;
		}
		double d1 = (d * d + r1 * r1 - r2 * r2) / 2.0d / d;
		double d2 = (d * d + r2 * r2 - r1 * r1) / 2.0d / d;
		s12[0] = 2.0d * Math.PI * r1 * (r1 - d1);
		s12[1] = 2.0d * Math.PI * r2 * (r2 - d2);
		return 0;
	} 
	
	public static final double calcSAVOL(double s) {
		return 0.877553507d * s + 17.64400201d;
	}
	
	public static final double calcConnoly(double s) {
		return 0.810491313d * s + 32.28203361d;
	}
	
	public static final double calcSpartan(double s) {
		return 0.919322875d * s + 17.10834414d;
	}
	
	public static Fragment[][] getHeteroRingsPSA(Molecule mol) throws SearchException, IOException {
		Rings.setMolecule(mol);
		int[][] heteroAliphaticRings = Rings.heteroAliphaticRings();
		int[][] heteroAromaticRings = Rings.heteroAromaticRings();
		Fragment[] heteroAliphaticRingsPSA = null;
		Fragment[] heteroAromaticRingsPSA = null;
		if(heteroAliphaticRings != null) {
			heteroAliphaticRingsPSA = new Fragment[heteroAliphaticRings.length];
			for(int i = 0; i < heteroAliphaticRings.length; i++) {
				heteroAliphaticRingsPSA[i] = new Fragment(mol, heteroAliphaticRings[i]);
				double[] psa = new double[2];
				getPSA(mol, psa, heteroAliphaticRings[i]);
				heteroAliphaticRingsPSA[i].addProp("PSA", new Double(psa[0]));
			}
			Arrays.sort(heteroAliphaticRingsPSA, new Comparator<Fragment>() {
				public int compare(Fragment f1, Fragment f2) {
					if(f1.getMolecule().getMass() > f2.getMolecule().getMass()) return 1;
					else if(f1.getMolecule().getMass() < f2.getMolecule().getMass()) return -1;
					else return 0;
				}
			});
		}
		if(heteroAromaticRings != null) {
			heteroAromaticRingsPSA = new Fragment[heteroAromaticRings.length];
			for(int i = 0; i < heteroAromaticRings.length; i++) {
				heteroAromaticRingsPSA[i] = new Fragment(mol, heteroAromaticRings[i]);
				double[] psa = new double[2];
				getPSA(mol, psa, heteroAromaticRings[i]);
				heteroAromaticRingsPSA[i].addProp("PSA", new Double(psa[0]));
			}
			Arrays.sort(heteroAromaticRingsPSA, new Comparator<Fragment>() {
				public int compare(Fragment f1, Fragment f2) {
					if(f1.getMolecule().getMass() > f2.getMolecule().getMass()) return 1;
					else if(f1.getMolecule().getMass() < f2.getMolecule().getMass()) return -1;
					else return 0;
				}
			});
		}
		Fragment[][] rings = new Fragment[2][];
		if(heteroAliphaticRingsPSA != null) {
			rings[0] = new Fragment[heteroAliphaticRingsPSA.length]; 
			for(int i = heteroAliphaticRingsPSA.length - 1; i >= 0; i--) {
				rings[0][rings[0].length - i - 1] = heteroAliphaticRingsPSA[i];
			}
		}
		if(heteroAromaticRingsPSA != null) {
			rings[1] = new Fragment[heteroAromaticRingsPSA.length]; 
			for(int i = heteroAromaticRingsPSA.length - 1; i >= 0; i--) {
				rings[1][rings[1].length - i - 1] = heteroAromaticRingsPSA[i];
			}
		}
		return rings;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws SearchException 
	 */
	public static void main(String[] args) throws IOException, SearchException {
		// TODO Auto-generated method stub
		Molecule mol = MolImporter.importMol("C1NC(C2CNC(C3CCNCC4=CC=NC=C34)C3=CNC=C23)C2=C/C=C\\N=C/C=C1\\2", "smiles");
		Hydrogenize.addHAtoms(mol);
		mol.aromatize(MoleculeGraph.AROM_GENERAL);
		Fragment[][] heteroRings = getHeteroRingsPSA(mol);
		// heteroaliphatic rings
		if(heteroRings[0] != null) {
			System.out.println("Ring\tPSA");
			for(Fragment hRing : heteroRings[0]) {
				System.out.println(hRing.toString("smiles") + "\t" + hRing.getProp("PSA"));
			}
		}
		// heteroaromatic rings
		if(heteroRings[1] != null) {
			System.out.println("Ring\tPSA");
			for(Fragment hRing : heteroRings[1]) {
				System.out.println(hRing.toString("smiles") + "\t" + hRing.getProp("PSA"));
			}
		}
	}

}
