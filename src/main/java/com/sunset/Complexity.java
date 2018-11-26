package com.sunset;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import chemaxon.core.util.BondTable;
import chemaxon.formats.MolImporter;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;

/**
 * Computes Synthetic and Molecular Complexity for in Silico Chemistry
 * according to the article in: J. Chem. Inf. Model., 45 (5), 1237 -1243, 2005. 10.1021/ci0501387
 * adapted from complexity.cpp based on OEChem toolkit by Tharun
 * @author Oleg Ursu
 *
 */
public class Complexity {
	
	public static final Pattern[] ringsPatterns = new Pattern[] {
		new Pattern("*~1~*~*~1"),
		new Pattern("*~1~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~*~*~*~*~1"),
		new Pattern("*~1~*~*~*~*~*~*~*~*~*~1")
	};
	
	public static final Pattern[] patterns2 = new Pattern[] {
		new Pattern("[CH2]-[CH2]"),
		new Pattern("[C!H0;C!H3]-[C!H0;C!H3]"),
		new Pattern("[$([#6]([F,Cl,Br,I])[F,Cl,Br,I])]"),
		new Pattern("[H]"),
		new Pattern("O=[C,N;R0]"),
		new Pattern("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"),
		new Pattern("[nH1]"),
		new Pattern("[N+,n+]"),
		new Pattern("[N-,n-]"),
		new Pattern("[O-,o-]"),
		new Pattern("[$([#6](-!@[#7])-!@[#7!H0])]"),
		new Pattern("[$([#6](@[#7])@[#7!H0])]"),
		new Pattern("[$([A!C!H1](-@[#6])-@[#6])]"),
		new Pattern("[$([#6](=[#8])-[#6]=[!#6])]"),
		new Pattern("[$([#16X3v4,#16X4v6]([#6])[A])]"),
		new Pattern("[$(C(~O)(~O)~O)]")
	};
	
	public static final Pattern[] patterns = new Pattern[] {
		new Pattern("[F]"),
		new Pattern("[Cl]"),
		new Pattern("[Br]"),
		new Pattern("[I]"),
		new Pattern("[H]"),
		new Pattern("[#5]"),
		new Pattern("[#6]"),
		new Pattern("[#16]"),
		new Pattern("[P]"),
		new Pattern("[$([C]-[C])]"),
		new Pattern("[$([C]=[C])]"),
		new Pattern("[$([C]#[C])]"),
		new Pattern("[$([C]-[N])]"),
		new Pattern("[#7]"),
		new Pattern("[$([N]=[N])]"),
		new Pattern("[$([C]-[O])]"),
		new Pattern("[$([C]-[F])]"),
		new Pattern("[$([C]-[P])]"),
		new Pattern("[$([C]-[S])]"),
		new Pattern("[#8]"),
		new Pattern("[$([C]-[Cl])]"),
		new Pattern("[$([C]-[Br])]"),
		new Pattern("[$([C]-[I])]"),
		new Pattern("[$([N]-[N])]"),
		new Pattern("[$([N]-[O])]"),
		new Pattern("[$([O]-[S])]"),
		new Pattern("[$([C!D1]-@C(=O)-@[N!a])]"),
		new Pattern("[$([C!D1]-!@C(=O)-!@[N!a])]"),
		new Pattern("[$([#6]=,:[*R!C!c]=,:[#6R]~@[*R!C!c])]"),
		new Pattern("[$([A!C!H1](-!@[#6])-!@[#6])]"),
		new Pattern("[$([#7]=[C!$(*a)]-!@[N!H0])]"),
		new Pattern("[$([#8]=[#6H0]-!@[N!H2])!$(NC[!#6])]"),
		new Pattern("[#8]=[#6](-!@[NH1,NH0])-!@[N!H2]"),
		new Pattern("[$(O=[CD3]([#6])[#6])]"),
		new Pattern("[$([#16X2v2]([#6])[#6])]"),
		new Pattern("[$([#8X2H0][#6]=[#8D1])]"),
		new Pattern("[$([#8X2v2]([#6])[#6])]"),
		new Pattern("c1:c-@C-@N-@C-@C-@O-@C-@c:c-1"),
		new Pattern("[N,O,C]C(=O)C1=C[C@H](*)C[C@H](O)O1"),
		new Pattern("C[C@H]1O[C@H]~3O[C@H]~C~2~C~C[C@H]([C@@H]~1)[C@@H]~23"),
		new Pattern("C[C@@H]1O[C@@H]~2O[C@H][C@@H][C@@H]~3~C~C~C~1[C@@H]~23"),
		new Pattern("C[C@H]2C[C@]14~C~C~C~C[C@H]1Oc3cccc(CN2)c34"),
		new Pattern("[#6][C@H]1C[C@@H]([#6])O[C@@H](-a)O1"),
		new Pattern("C2=CC[C@@H]1C(=O)~*~*~C(=O)[C@@H]1[C@@H]2-a"),
		new Pattern("C2=CC[C@@H]1C(=O)~*~C(=O)[C@@H]1[C@@H]2-a"),
		new Pattern("a-[CH,CH2;R0;0*]"),
		new Pattern("[R;0*]-[CH2R0,NHR0,OR0;0*]-[R]"),
		new Pattern("*-[CD3H,ND2;R0;0*](-a)-a"),
		new Pattern("[a]-&!@[a]"),
		new Pattern("[NR;0*]-[CD3R0;0*](=O)-[R]"),
		new Pattern("[NR;0*]-[CD2R0;0*]-[R]"),
		new Pattern("[NR;0*]-[CD2R0;0*]-[CD2,CD3,OD2,ND2,ND3,aD2,aD3]"),
		new Pattern("a-[NHR0]-[CR0;0*](=O)-[OR0,NR0,0*]"),
		new Pattern("[CR,NR]=[CR]-&!@[a]"),
		new Pattern("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"),
		new Pattern("O=[C,N;R0]"),
		new Pattern("[O,N;!H0]"),
		new Pattern("[C]([R])([R])([R])[R]"),
		new Pattern("[$([#6X4@](*)(*)(*)*),$([#6X4@H](*)(*)*)]"),
		new Pattern("[$([#6X4@](*)(*)(*)*),$([#6X4@H](*)(*)*)][$([#6X4@](*)(*)(*)*),$([#6X4@H](*)(*)*)]"),
		new Pattern("[$([#6]([F,Cl,Br,I])[F,Cl,Br,I])]"),
		new Pattern("[C!H0;C!H3]-[C!H0;C!H3]"),
		new Pattern("[CH2]-[CH2]"),
		new Pattern("[$([C]=[N])]"),
		new Pattern("[$([C]=[O])]"),
		new Pattern("[$([C]=[S])]"),
		new Pattern("[$([N]=[O])]"),
		new Pattern("[$([O]=[S])]"),
		new Pattern("[$([C]#[N])]"),
		new Pattern("[$([#6]:[#6])]"),
		new Pattern("[$([#6]:[#7])]"),
		new Pattern("[$([#6]:[#16])]"),
		new Pattern("[$([#7]:[#7])]"),
		new Pattern("[$([#7]:[#8])]")
	};
	
	private static final int[] metric = new int[90];
	private static final int[] metric2 = new int[26];
	
	private static int[][] ringConnections;
	private static int[][] rings;
	private static int[][] ringSystems;
	private static int[] ringSystemsAtomIdx;
	
	private static int bridged = 0;
	private static int ringComplexity = 0;
	private static int flag = 0;
	
	private static final void reset() {
		for(int i = 0; i < metric.length; i++) {
			metric[i] = 0;
		}
		for(int i = 0; i < metric2.length; i++) {
			metric2[i] = 0;
		}
		ringConnections = null;
		rings = null;
		ringSystems = null;
		ringSystemsAtomIdx = null;
	}
	
	private static final void init(Molecule mol) {
		rings = mol.getSSSR();
		findRingConnections();
		findRingSystems();
		assignRingSystemAtomIdx(mol);
	}
	
	private static final int rings(Molecule mol) throws SearchException {
		bridged = 0;
		ArrayList<Integer> ringSizes = new ArrayList<Integer>();
		int rct = -1;
		ArrayList<int[]> ringBonds = new ArrayList<int[]>();
		for(Pattern patt : ringsPatterns) {
			patt.pattern.setTarget(mol);
			int[][] hits = patt.pattern.findAll();
			if(hits != null) {
				for(int[] hit : hits) {
					ringSizes.add(hit.length);
					ringBonds.add(hitBonds(hit, mol, patt.pattern.getQuery()));
					rct++;
				}
			}
		}
		int samebond = 0;
		int bridgedct = 0;
		int totbridged = 0;
		for(int ii = 0; ii <= rct; ii++) {
			for(int jj = 0; jj<= rct && jj != ii; jj++) {
				samebond = 0;
				for(int kk = 0; kk < ringBonds.get(ii).length; kk++) {
					for(int ll = 0; ll < ringBonds.get(jj).length; ll++) {
						if(ringBonds.get(ii)[kk] == ringBonds.get(jj)[ll]) {
							samebond++;
						}
					}
				}
				if(samebond >= 2 && samebond != ringSizes.get(ii) - 1 && samebond != ringSizes.get(jj) - 1) {
					bridgedct++;
				}
			}
			if(bridgedct > 0) {
				totbridged++;
			}
		}

		totbridged /= 2;
		return totbridged;
	}
	
	public static final int[] hitBonds(int[] hit, Molecule target, Molecule query) {
		int[] bonds = new int[hit.length];
		int n = query.getBondCount();
		BondTable bondTable = target.getBondTable();
		for(int i = 0; i < n; i++) {
			MolBond bond = query.getBond(i);
			int qid1 = query.indexOf(bond.getAtom1());
			int qid2 = query.indexOf(bond.getAtom2());
			bonds[i] = bondTable.getBondIndex(hit[qid1], hit[qid2]);
		}
		return bonds;
	}

	
	public static final double complexity(Molecule mol) throws SearchException {
		reset();
		init(mol);
		int z = 0;
		int[][] matches = null;
		for(int i = 0; i < patterns.length; i++) {
			patterns[i].pattern.setTarget(mol);
			matches = patterns[i].pattern.findAll();
			if(matches != null) {
				metric[z] = matches.length;
			}
			z++;
		}
		
		int z2 = 0;
		for(int i = 0; i < patterns2.length; i++) {
			patterns2[i].pattern.setTarget(mol);
			matches = patterns2[i].pattern.findAll();
			if(matches != null) {
				metric2[z2] = matches.length;
			}
			z2++;
		}
		
		int hCount = 0;
		
		for(int i = 0; i < mol.getAtomCount(); i++) {
			MolAtom atom = mol.getAtom(i);
			if(atom.getAtno() == 1) {
				hCount++;
			} else {
				hCount += atom.getExplicitHcount();
				hCount += atom.getImplicitHcount();
			}
		}
		
		int ringCt = ringSystems.length;
		int[] atomRingSize = smallestAtomRingSize(mol);
		int three=0,four=0,five=0,six=0,seven=0,eight=0,nine=0;
		for(int ringSize : atomRingSize) {
			switch (ringSize) {
			case 3:
				three++;
				break;
			case 4:
				four++;
				break;
			case 5:
				five++;
				break;
			case 6:
				six++;
				break;
			case 7:
				seven++;
				break;
			case 8:
				eight++;
				break;
			case 9:
				nine++;
				break;
			}
		}
		three = (int)Math.ceil(three/3.0d);
		four = (int)Math.ceil(four/4.0d);
		five = (int)Math.ceil(five/5.0d);
		six = (int)Math.ceil(six/6.0d);
		seven = (int)Math.ceil(seven/7.0d);
		eight = (int)Math.ceil(eight/8.0d);
		nine = (int)Math.ceil(nine/9.0d);
		int nor = mol.getBondCount() - mol.getAtomCount() + 1;
		metric[z++] = nor;
		ringComplexity=(three+four+seven+eight+nine)*2+(five+six)*1;
		int totringSystems = ringCt;
		metric[z++] = totringSystems;
		if(nor > 2 && nor <= 10) {
			bridged = countBridged(mol);
		} else {
			bridged = 0;
			if(nor > 10) {
				System.err.println("WARRNING: This version cannot find bridges for molecules with more than 10 rings.");
			}
		}
		int tfused = totfused(mol);
		metric[z++] = tfused;
		metric[4] = hCount;
		metric[6] -= (metric[10] + metric[11]);
		HashSet<Integer> adjChiralIndex = new HashSet<Integer>();
		int perceivedChiral = 0;
		for(int i = 0; i < mol.getAtomCount(); i++) {
			if(mol.getParity(i) != 0) {
				perceivedChiral++;
				adjChiralIndex.add(i);
			}
		}
		if(metric[58] < perceivedChiral) {
			metric[58] = perceivedChiral;
			flag = 1;
		}
		int adjChiralPerceived=0, adjChiralCount=0;
		int[][] btab = mol.getBondTable().getMatrixArray();
		for(Integer index : adjChiralIndex) {
			for(int j = 0; j < btab[index].length; j++) {
				if(btab[index][j] != -1 && j != index && adjChiralIndex.contains(j)) {
					adjChiralCount = 1;
				}
			}
			if(adjChiralCount == 1) {
				adjChiralCount = 0;
				adjChiralPerceived++;
			}
		}
		if(metric[59] < adjChiralPerceived) {
			metric[59] = adjChiralPerceived;
		}
		if((metric[58]>metric[59]) && metric[58]!=1 && metric[59]!=0) {
			metric[59] += 1;
		}
		return computeComplexity(metric, metric2);
	}
	
	public static final double computeComplexity(int[] m, int[] m2) {
		double comp = 0D;
		
		comp=m[0] * 1.446 +
		m[1] * 1.384 + 
		m[2] * 1.244 + 
		m[3] *1.103 + 
        m[4] * 0.0 + 
        m[5] * 0.851 + 
        m[6] * 1.0 + 
        m[7] * 1.235 + 
        m[8] * 1.086+ 
		m[9] * 1.0 + 
		m[12] * 0.857 +
		m[10] * 0.5 +
		m[11] * 0.333 +
        m[13] * 1.149 + 
        m[14] * 0.367 +
		m[15] * 0.750 +
        m[16] * 0.667 +
		m[17] * 0.400 +
		m[18] * 0.375 +
        m[19] * 1.297 +
        m[20] * 0.353 +
		m[21] * 0.171 +
		m[22]  * 0.113 +
		m[23] * 0.735 +
		m[24] * 0.643 +
		m[25] * 0.281 -
        2.0* m[26]-2.0* m[27]-2.0* m[28]-2.0* m[29]-2.0* m[30]-2.0* m[31]-2.0* m[32]-2.0* m[33]-2.0* m[34]-2.0* m[35] -
        2.0* m[36]-2.0* m[37]-2.0* m[38]-2.0* m[39]-2.0* m[40]-2.0* m[41]-2.0* m[42]-2.0* m[43]-2.0* m[44]-2.0* m[45] -
        2.0* m[46]-2.0* m[47]-2.0* m[48]-2.0* m[49]-2.0* m[50]-2.0* m[51]-2.0* m[52]-2.0*m [53] -
        m[54] * 0.0 +
        m[56] * 0.0 +
        m[57] * 3.0 +
        m[58] * 2.0 -
        m[59] * 1.0 +
        m[60] * 0.0 +
        m[75] * 1 -
		m[61] * 0.0 -
		m[62] * 0.0 +
		m[63] * 0.429 +
		m[64] * 0.375 +
		m[65] * 0.188 +
		m[66] * 0.321 +
		m[67] * 0.141 +
		m[68] * 0.286 +
		m[69] * 0.667 +
		m[70] * 0.571 +
		m[71] * 0.250 +
		m[72] * 0.490 +
		m[73] * 0.423 +
		ringComplexity;	
		
		if(bridged > 0) {
			comp += 1*4.0D;
		}
		
		comp += m[76]*2.0D;
		if(flag == 1) {
			flag = 0;
			comp -= 1;
		}
		
		comp += m2[0] * 1.0 * 0.0 * +
        		m2[1] * 1.0 * 0.0 +
        		m2[2] * 3.0 * 0.0 +
        		m2[3] * 0.0 * 0.0 + 
        		m2[4] * 1.0 * 0.0 +
        		m2[5] * 1.0 * 0.0 +
        		m2[6] * 2.0 * 0.0 +
        		(m2[7] + m2[8] + m2[9] ) * 3.0 * 0.0 -
        		2.0 *( m2[10] + m2[11] + m2[12] + m2[13] + m2[14] + m2[15]);
		return comp;
	}
	
	public static final int countBridged(Molecule mol) throws SearchException {
		return rings(mol);
	}
	
	public static final int[] smallestAtomRingSize(Molecule mol) {
		int[] atomRingSize = new int[mol.getAtomCount()];
		for(int[] ring : rings) {
			for(int atom : ring) {
				if(atomRingSize[atom] == 0 || atomRingSize[atom] > ring.length) {
					atomRingSize[atom] = ring.length;
				}
			}
		}
		return atomRingSize;
	}
	
	private static final void findRingConnections() {
		ringConnections = new int[rings.length][rings.length];
        for(int i = 0; i < rings.length; i++) {
            for(int j = 0; j < rings[i].length; j++) {
                for(int m = i + 1; m < rings.length; m++) {
                    for(int n = 0; n < rings[m].length; n++)
                    	if(rings[i][j] == rings[m][n]) {
                            ringConnections[i][m]++;
                            ringConnections[m][i]++;
                        }
                }
            }
        }
	}
	
	private static final void findRingSystems() {
		int systems[] = new int[rings.length];
        int c = 0;
        for(int i = 0; i < ringConnections.length; i++) {
            if(systems[i] == 0)
                systems[i] = ++c;
            for(int j = i + 1; j < ringConnections[i].length; j++) {
                if(ringConnections[i][j] <= 0)
                	continue;
                if(systems[j] == 0) {
                    systems[j] = systems[i];
                    continue;
                }
                int m = systems[i];
                for(int k = 0; k < systems.length; k++) {
                    if(systems[k] == m) {
                        systems[k] = systems[j];
                    }
                }
            }
        }

        boolean r[] = new boolean[systems.length];
        for(int i = 0; i < systems.length; i++) {
            r[systems[i] - 1] = true;
        }
        
        c = 0;
        for(int i = 0; i < r.length; i++) {
            if(r[i])
                c++;
        }
        ringSystems = new int[c][];
        for(int i = 0; i < c; i++) {
            int d = 0;
            for(int j = 0; j < systems.length; j++)
                if(systems[j] == i + 1)
                    d++;

            ringSystems[i] = new int[d];
            d = 0;
            for(int j = 0; j < systems.length; j++)
                if(systems[j] == i + 1) {
                    ringSystems[i][d] = j;
                    d++;
                }
        }
	}
	
	
	/*
	 * assign ring system id to each atom in molecule
	 */
	
	private static final void assignRingSystemAtomIdx(Molecule mol) {
		ringSystemsAtomIdx = new int[mol.getAtomCount()];
		for(int i = 0; i < ringSystems.length; i++) {
			for(int ring : ringSystems[i]) {
				for(int idx : rings[ring]) {
					ringSystemsAtomIdx[idx] = i + 1;
				}
			}
		}
	}
	
	private static final int totfused(Molecule mol) {
		if(rings.length == 0) {
			return 0;
		}
		int[] atomCt = new int[ringSystems.length + 1];
		int[] bondCt = new int[ringSystems.length + 1];
		int[] spiroCt = new int[ringSystems.length + 1];
		int ringIdx, ringDegree, idx;
		int totct = 0;
		metric[57] = 0;
		int[][] btab =  mol.getBondTable().getMatrixArray();
		for(idx = 0; idx < mol.getAtomCount(); idx++) {
			ringIdx = ringSystemsAtomIdx[idx];
			atomCt[ringIdx]++;
			ringDegree = 0;
			for(int bondIdx : btab[idx]) {
				if(bondIdx != -1 && mol.isRingBond(bondIdx)) {
					ringDegree++;
				}
			}
			if(ringDegree > 3) {
				spiroCt[ringIdx] += ringDegree - 3;
			}
		}
		for(int bond = 0; bond < mol.getBondCount(); bond++) {
			idx = mol.indexOf(mol.getBond(bond).getAtom1());
			ringIdx = ringSystemsAtomIdx[idx];
			if(ringIdx > 0 && ringIdx == ringSystemsAtomIdx[mol.indexOf(mol.getBond(bond).getAtom2())]) {
				bondCt[ringIdx]++;
			}
		}
		for(int i = 1; i <= ringSystems.length; i++) {
			int ct = bondCt[i] - atomCt[i] + 1 - spiroCt[i];
			if(ct > 1) {
				totct += ct;
			}
		}
		metric[57] = 0;
		for(int i = 0; i <= ringSystems.length; i++) {
			metric[57] += spiroCt[i];
		}
		return totct;
	}
	
	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException, SearchException {
		Molecule mol = MolImporter.importMol("C[C@@](O)(CCO)CC(O)=O", "smiles");
		mol.aromatize(MoleculeGraph.AROM_GENERAL);
		System.out.println(complexity(mol));
	}

}
