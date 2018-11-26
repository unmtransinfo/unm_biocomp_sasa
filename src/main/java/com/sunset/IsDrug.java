package com.sunset;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.kohsuke.args4j.ParserProperties;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.io.MolExportException;
import chemaxon.reaction.Standardizer;
import chemaxon.reaction.StandardizerException;
import chemaxon.sss.SearchConstants;
import chemaxon.sss.search.MolSearchOptions;
import chemaxon.sss.search.SearchException;
import chemaxon.sss.search.StandardizedMolSearch;
import chemaxon.struc.Molecule;

/**
 * @author Oleg Ursu
 * label drugs in input file
 */
public class IsDrug {
	
	private static final Logger logger = LoggerFactory.getLogger(IsDrug.class);

	@Argument(metaVar = "[drugnames.txt] [drugdb.sdf] [input] [output]", required = true, usage = "drugnames/drugdb/input/output file name, default none")
	private List<String> arguments = new ArrayList<>();
	
	@Option(name = "-h", aliases = {"--help"}, required = false, usage = "print this message")
	private boolean help = false;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		IsDrug app = new IsDrug();
		ParserProperties properties = ParserProperties.defaults();
        properties.withUsageWidth(200);
        CmdLineParser parser = new CmdLineParser(app, properties);
		try {
			parser.parseArgument(args);
			if(app.help) {
				parser.printUsage(System.out);
				return;
			}
			if(app.arguments.size() != 4) {
				parser.printUsage(System.out);
				return;
			}
		} catch (CmdLineException e) {
			// TODO Auto-generated catch block
			logger.error(e.getMessage());
		}
		app.run();
	}
	
	public Set<String> readDrugNames() {
		Set<String> set = new HashSet<>();
		try(BufferedReader reader = Files.newBufferedReader(Paths.get(arguments.get(0)), Charset.defaultCharset())) {
			String line;
			while((line = reader.readLine()) != null) {
				set.add(line.trim().toLowerCase());
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(e.getMessage(), e);
		}
		return set;
	}
	
	@SuppressWarnings("deprecation")
	public List<Molecule> readDrugList() {
		MolImporter reader = null;
		List<Molecule> list = null;
		Molecule mol = null;
		try {
			reader = new MolImporter(arguments.get(1));	
			list = new ArrayList<>();
			while((mol = reader.read()) != null) {
				list.add(mol);
			}
		} catch (MolFormatException e) {
			// TODO Auto-generated catch block
			logger.error(e.getMessage());
			logger.error("last good mol id {}", mol.getProperty("ID"));
			return null;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(e.getMessage());
			return null;
		} finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					logger.error(e.getMessage());
				}
			}
		}
		if(list.size() > 0) {
			return list;
		} else {
			return null;
		}
	}
	
	@SuppressWarnings("deprecation")
	public void run() {
		List<Molecule> drugList = readDrugList();
		if(drugList == null || drugList.size() == 0) {
			logger.error("reading drug list from {}", arguments.get(0));
			return;
		}
		MolImporter reader = null;
		MolExporter writer = null;
		try {
			 reader = new MolImporter(arguments.get(2));
			 writer = new MolExporter(arguments.get(3), "rdf");
			 MolSearchOptions mso = new MolSearchOptions(SearchConstants.DUPLICATE);
			 mso.setStereoSearchType(SearchConstants.STEREO_IGNORE);
			 StandardizedMolSearch ms = new StandardizedMolSearch();
			 Standardizer std = new Standardizer(this.getClass().getResourceAsStream("/standardizer.xml"));
			 ms.setSearchOptions(mso);
			 ms.setStandardizer(std, true, true);
			 Molecule mol;
			 boolean match;
			 Set<String> names = readDrugNames();
			 while((mol = reader.read()) != null) {
				 if(mol.getProperty("ROOT:MOL.GNAME") != null && names.contains(mol.getProperty("ROOT:MOL.GNAME").toLowerCase())) {
					 mol.setProperty("ROOT:MOL.ISDRUG", "1");
					 writer.write(mol);
					 continue;
				 }
				 if(mol.getProperty("ROOT:MOL.NAME") != null && names.contains(mol.getProperty("ROOT:MOL.NAME").toLowerCase())) {
					 mol.setProperty("ROOT:MOL.ISDRUG", "1");
					 writer.write(mol);
					 continue;
				 }
				 match = false;
				 ms.setTarget(mol);
				 for(Molecule query : drugList) {
					 if(!query.getFormula().equals(mol.getFormula())) {
						 continue;
					 }
					 ms.setQuery(query);
					 try {
						if(ms.isMatching()) {
							 match = true;
							 break;
						 }
					} catch (SearchException e) {
						// TODO Auto-generated catch block
						logger.error(e.getLocalizedMessage());
					}
				 }
				 if(match) {
					 mol.setProperty("ROOT:MOL.ISDRUG", "1");
				 } else {
					 mol.setProperty("ROOT:MOL.ISDRUG", "0");
				 }
				 writer.write(mol);
			 }
		} catch (MolFormatException e) {
			// TODO Auto-generated catch block
			logger.error(e.getMessage());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error(e.getMessage());
		} catch (StandardizerException e) {
			// TODO Auto-generated catch block
			logger.error(e.getLocalizedMessage());
		} finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					logger.error(e.getMessage());
				}
			}
			if(writer != null) {
				try {
					writer.close();
				} catch (MolExportException e) {
					// TODO Auto-generated catch block
					logger.error(e.getMessage());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					logger.error(e.getMessage());
				}
			}
		}
		
	}
}
