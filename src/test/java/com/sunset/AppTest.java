package com.sunset;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;


/**
 * Unit test for simple App.
 */
public class AppTest 
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    /**
     * Rigourous Tests :-)
     */
    private static final String testsmi = "NCCc1ccc(O)c(O)c1";
    public void testComplexity() throws Exception
    {
	Molecule mol = MolImporter.importMol(testsmi, "smiles");
	mol.aromatize(MoleculeGraph.AROM_GENERAL);
	double cmplx = Complexity.complexity(mol);
        assertTrue( cmplx > 0.0 );
    }
    public void testABE() throws Exception
    {
	Molecule mol = MolImporter.importMol(testsmi, "smiles");
	mol.aromatize(MoleculeGraph.AROM_GENERAL);
	double abe = ABE.calcABE(mol);
        assertTrue( abe > 0.0 );
    }
    public void testPSA() throws Exception
    {
	Molecule mol = MolImporter.importMol(testsmi, "smiles");
	mol.aromatize(MoleculeGraph.AROM_GENERAL);
	double[] psa = new double[2];
	PSA.getPSA(mol, psa);
        double sfa_perc_pol=psa[0]/(psa[0]+psa[1])*100.0;
        double sfa_perc_nonpol=psa[1]/(psa[0]+psa[1])*100.0;
        assertTrue( (sfa_perc_pol>0.0) && (sfa_perc_nonpol>0.0) );
    }
}
