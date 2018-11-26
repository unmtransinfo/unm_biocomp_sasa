package com.sunset;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Types;
import java.util.HashMap;
import java.util.Properties;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.io.MPropHandler;
import chemaxon.struc.Molecule;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

/**
 * @author Oleg Ursu
 * Import SDF file into ORACLE database
 */
public class OracleImportSDF {

	@Parameter(names = {"-m", "--mapping"}, description = "mapping file", required = true)
	private String mappingFileName;
	
	@Parameter(names = {"-c", "--config"}, description = "configuration file", required = true)
	private String configFileName;
	
	@Parameter(names = {"-i", "--input"}, description = "input file", required = true)
	private String inputFileName;
		
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		OracleImportSDF app = new OracleImportSDF();
		new JCommander(app, args);
		app.run();
	}
	
	public void run() {
		Properties config = getConfiguration();
		if(config == null || config.size() == 0) {
			System.out.println("Error reading configuration");
			return;
		}
		Connection conn = getConnection(config);
		if(conn == null) {
			System.out.println("Falied to establish connection to " + config.getProperty("jdbcUrl"));
			return;
		}
		HashMap<String, String> mapping = getMapping();
		if(mapping == null || mapping.size() == 0) {
			System.out.println("Error reading mapping");
			return;
		}
		String insertStmt = "INSERT INTO " + mapping.get("__TABLE_NAME__") + "(";
		insertStmt += mapping.get("SMDL.ID");
		insertStmt += "," + mapping.get("SMDL.IDx");
		insertStmt += "," + mapping.get("SMDL.SID");
		insertStmt += "," + mapping.get("MOL.REF");
		insertStmt += "," + mapping.get("MOL.NAME");
		insertStmt += "," + mapping.get("MOL.GNAME");
		insertStmt += "," + mapping.get("MOL.SMI");
		insertStmt += "," + mapping.get("EST.LOGKOW");
		insertStmt += "," + mapping.get("EST.LOGWSOL");
		insertStmt += "," + mapping.get("NO.ATOMS");
		insertStmt += "," + mapping.get("NO.BONDS");
		insertStmt += "," + mapping.get("NO.RINGS");
		insertStmt += "," + mapping.get("NO.ROT.BONDS");
		insertStmt += "," + mapping.get("NO.RIG.BONDS");
		insertStmt += "," + mapping.get("NO.HETERO.ATOMS");
		insertStmt += "," + mapping.get("NO.NONPOL.ATOMS");
		insertStmt += "," + mapping.get("NO.POS.IONIZ");
		insertStmt += "," + mapping.get("NO.NEG.IONIZ");
		insertStmt += "," + mapping.get("LPK.HB.DON");
		insertStmt += "," + mapping.get("LPK.HB.ACC");
		insertStmt += "," + mapping.get("LPK.CLOGP");
		insertStmt += "," + mapping.get("LPK.SCORE");
		insertStmt += "," + mapping.get("SFA.POL");
		insertStmt += "," + mapping.get("SFA.NONPOL");
		insertStmt += "," + mapping.get("SFA.PERC.POL");
		insertStmt += "," + mapping.get("SFA.PERC.NONPOL");
		insertStmt += "," + mapping.get("MOL.XMR");
		insertStmt += "," + mapping.get("MOL.ABE");
		insertStmt += "," + mapping.get("MOL.SMCM");
		insertStmt += "," + mapping.get("MOL.ISDRUG");
		insertStmt += "," + mapping.get("NO.AROMATIC.RINGS");
		insertStmt += "," + mapping.get("NO.ALIPHATIC.RINGS");
		insertStmt += "," + mapping.get("__STRUCTURE_COLUMN__");
		insertStmt += ") VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
		try {
			PreparedStatement ps = conn.prepareStatement(insertStmt);
			MolImporter molReader = new MolImporter(inputFileName);
			Molecule mol;
			long count = 0L;
			float f = 0.0F;
			int i = 0;
			short s = 0;
			String str = null;
			while((mol = molReader.read()) != null) {
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "SMDL.ID"));
					ps.setInt(1, i);
				} catch(Exception e) {
					ps.setNull(1, Types.INTEGER);
				}
				str = MPropHandler.convertToString(mol.properties(), "SMDL.IDx");
				if(str == null || str.length() == 0) {
					ps.setNull(2, Types.CHAR);
				} else {
					ps.setString(2, str);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "SMDL.SID"));
					ps.setInt(3, i);
				} catch(Exception e) {
					ps.setNull(3, Types.INTEGER);
				}
				str = MPropHandler.convertToString(mol.properties(), "MOL.REF");
				if(str == null || str.length() == 0) {
					ps.setNull(4, Types.VARCHAR);
				} else {
					ps.setString(4, str);
				}
				str = MPropHandler.convertToString(mol.properties(), "MOL.NAME");
				if(str == null || str.length() == 0) {
					ps.setNull(5, Types.VARCHAR);
				} else {
					ps.setString(5, str);
				}
				str = MPropHandler.convertToString(mol.properties(), "MOL.GNAME");
				if(str == null || str.length() == 0) {
					ps.setNull(6, Types.VARCHAR);
				} else {
					ps.setString(6, str);
				}
				str = MPropHandler.convertToString(mol.properties(), "MOL.SMI");
				if(str == null || str.length() == 0) {
					ps.setNull(7, Types.VARCHAR);
				} else {
					ps.setString(7, str);
				}
				try {
					f = Float.parseFloat(MPropHandler.convertToString(mol.properties(), "EST.LOGKOW"));
					ps.setFloat(8, f);
				} catch(Exception e) {
					System.out.println("Error parsing float from " + MPropHandler.convertToString(mol.properties(), "EST.LOGKOW"));
					ps.setNull(8, Types.FLOAT);
				}
				try {
					f = Float.parseFloat(MPropHandler.convertToString(mol.properties(), "EST.LOGWSOL"));
					ps.setFloat(9, f);
				} catch(Exception e) {
					ps.setNull(9, Types.FLOAT);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.ATOMS"));
					ps.setInt(10, i);
				} catch(Exception e) {
					ps.setNull(10, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.BONDS"));
					ps.setInt(11, i);
				} catch(Exception e) {
					ps.setNull(11, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.RINGS"));
					ps.setInt(12, i);
				} catch(Exception e) {
					ps.setNull(12, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.ROT.BONDS"));
					ps.setInt(13, i);
				} catch(Exception e) {
					ps.setNull(13, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.RIG.BONDS"));
					ps.setInt(14, i);
				} catch(Exception e) {
					ps.setNull(14, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.HETERO.ATOMS"));
					ps.setInt(15, i);
				} catch(Exception e) {
					ps.setNull(15, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.NONPOL.ATOMS"));
					ps.setInt(16, i);
				} catch(Exception e) {
					ps.setNull(16, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.POS.IONIZ"));
					ps.setInt(17, i);
				} catch(Exception e) {
					ps.setNull(17, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.NEG.IONIZ"));
					ps.setInt(18, i);
				} catch(Exception e) {
					ps.setNull(18, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "LPK.HB.DON"));
					ps.setInt(19, i);
				} catch(Exception e) {
					ps.setNull(19, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "LPK.HB.ACC"));
					ps.setInt(20, i);
				} catch(Exception e) {
					ps.setNull(20, Types.INTEGER);
				}
				try {
					f = Float.parseFloat(MPropHandler.convertToString(mol.properties(), "LPK.CLOGP"));
					ps.setFloat(21, f);
				} catch(Exception e) {
					ps.setNull(21, Types.FLOAT);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "LPK.SCORE"));
					ps.setInt(22, i);
				} catch(Exception e) {
					ps.setNull(22, Types.INTEGER);
				}
				try {
					f = Float.parseFloat(MPropHandler.convertToString(mol.properties(), "SFA.POL"));
					ps.setFloat(23, f);
				} catch(Exception e) {
					ps.setNull(23, Types.FLOAT);
				}
				try {
					f = Float.parseFloat(MPropHandler.convertToString(mol.properties(), "SFA.NONPOL"));
					ps.setFloat(24, f);
				} catch(Exception e) {
					ps.setNull(24, Types.FLOAT);
				}
				try {
					f = Float.parseFloat(MPropHandler.convertToString(mol.properties(), "SFA.PERC.POL"));
					ps.setFloat(25, f);
				} catch(Exception e) {
					ps.setNull(25, Types.FLOAT);
				}
				try {
					f = Float.parseFloat(MPropHandler.convertToString(mol.properties(), "SFA.PERC.NONPOL"));
					ps.setFloat(26, f);
				} catch(Exception e) {
					ps.setNull(26, Types.FLOAT);
				}
				try {
					f = Float.parseFloat(MPropHandler.convertToString(mol.properties(), "MOL.XMR"));
					ps.setFloat(27, f);
				} catch(Exception e) {
					ps.setNull(27, Types.FLOAT);
				}
				try {
					f = Float.parseFloat(MPropHandler.convertToString(mol.properties(), "MOL.ABE"));
					ps.setFloat(28, f);
				} catch(Exception e) {
					ps.setNull(28, Types.FLOAT);
				}
				try {
					f = Float.parseFloat(MPropHandler.convertToString(mol.properties(), "MOL.SMCM"));
					ps.setFloat(29, f);
				} catch(Exception e) {
					ps.setNull(29, Types.FLOAT);
				}
				try {
					s = Short.parseShort(MPropHandler.convertToString(mol.properties(), "MOL.ISDRUG"));
					ps.setShort(30, s);
				} catch(Exception e) {
					ps.setNull(30, Types.SMALLINT);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.AROMATIC.RINGS"));
					ps.setInt(31, i);
				} catch(Exception e) {
					ps.setNull(31, Types.INTEGER);
				}
				try {
					i = Integer.parseInt(MPropHandler.convertToString(mol.properties(), "NO.ALIPHATIC.RINGS"));
					ps.setInt(32, i);
				} catch(Exception e) {
					ps.setNull(32, Types.INTEGER);
				}
				str = MolExporter.exportToFormat(mol, "mol");
				if(str == null || str.length() == 0) {
					ps.setNull(33, Types.CLOB);
				} else {
					ps.setString(33, str);
				}
				ps.executeUpdate();
				count++;
				if(count % 100 == 0) {
					System.out.println(count + " records imported");
				}
			}
			ps.close();
			molReader.close();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			try {
				conn.commit();
				conn.close();
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
	}
	
	private Connection getConnection(Properties props) {
		Connection conn;
		try {
			conn = DriverManager.getConnection(props.getProperty("jdbcUrl"), props.getProperty("jdbcUser"), props.getProperty("jdbcPassword"));
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		return conn;
	}
	
	private Properties getConfiguration() {
		Properties props = new Properties();
		try {
			props.load(new FileInputStream(configFileName));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
		return props;
	}
	
	private HashMap<String, String> getMapping() {
		HashMap<String, String> map = new HashMap<String, String>();
		BufferedReader reader = null;
		String line;
		try {
			reader = new BufferedReader(new FileReader(mappingFileName));
			while((line = reader.readLine()) != null) {
				String[] tokens = line.split("=");
				map.put(tokens[0], tokens[1]);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		} finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
		}
		return map;
	}

}
