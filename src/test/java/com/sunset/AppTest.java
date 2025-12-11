package com.sunset;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;
import chemaxon.struc.MoleculeGraph;


/**
 * Unit test for simple App.  WITHOUT LICENSE MUST BE DISABLED.
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
  public void testDummy() throws Exception
  {
    assertTrue(true);
  }
  /********************************************************************
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
  public void testRingAliPSA() throws Exception
  {
    Molecule mol = MolImporter.importMol("CCC1=CN=C(C=C1)CCOC2=CC=C(C=C2)CC3C(=O)NC(=O)S3", "smiles");
    mol.aromatize(MoleculeGraph.AROM_GENERAL);
    Fragment[][] hetRings=null;
    hetRings=PSA.getHeteroRingsPSA(mol);
    double rpsa = 0.0;
    if (hetRings[0]!=null)
    {
    	rpsa = new Double(hetRings[0][0].getProp("PSA").toString());
    }
    assertTrue( rpsa > 0.0 );
  }
  public void testRingAroPSA() throws Exception
  {
    Molecule mol = MolImporter.importMol("CCC1=CN=C(C=C1)CCOC2=CC=C(C=C2)CC3C(=O)NC(=O)S3", "smiles");
    mol.aromatize(MoleculeGraph.AROM_GENERAL);
    Fragment[][] hetRings=null;
    hetRings=PSA.getHeteroRingsPSA(mol);
    double rpsa = 0.0;
    if (hetRings[1]!=null)
    {
    	rpsa = new Double(hetRings[1][0].getProp("PSA").toString());
    }
    assertTrue( rpsa > 0.0 );
  }
  *******************************************************************/
}
