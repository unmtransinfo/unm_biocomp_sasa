package edu.unm.health.biocomp.uniprot;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
/*

import uk.ac.ebi.kraken.interfaces.uniprot.Keyword;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.comments.Comment;
import uk.ac.ebi.kraken.interfaces.uniprot.comments.CommentType;
import uk.ac.ebi.kraken.interfaces.uniprot.description.Name;
import uk.ac.ebi.kraken.uuw.services.remoting.EntryRetrievalService;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtJAPI;
*/
/**
 * @author Oleg Ursu
 * Query UniProt database using UniProtJAPI {@link} http://www.ebi.ac.uk/uniprot/remotingAPI/
 */
public class UniProtQueryUtil {
	/*
	static {
		System.setProperty("org.apache.commons.logging.Log", "org.apache.commons.logging.impl.NoOpLog");
	}
		
	*/
	/*	

	private EntryRetrievalService entryRetrievalService = UniProtJAPI.factory.getEntryRetrievalService();
	private UniProtEntry entry = null;
	*/
	/**
	 * Query UniProtKB using primary accession numbers
	 * @param acc accession number as string
	 * @return true if an entry was found, false otherwise
	 */
		/*
	public boolean queryByAccessionNumber(String acc) {
		entry = entryRetrievalService.getUniProtEntry(acc);
		return entry != null;
	}
	*/
	
	/**
	 * Get scientific organism name for entry queried with one of the query methods 
	 * @return scientific organism name
	 */
			/*
	public String getOrganismScientificName() {
		return entry.getOrganism().getScientificName().toString();
	}
	
	public String getNCBITaxId() {
		return entry.getNcbiTaxonomyIds().get(0).getValue();
	}
		
	public String getProteinRecommendedName() {
		if(entry.getProteinDescription().hasRecommendedName()) {
			return entry.getProteinDescription().getRecommendedName().getFields().get(0).getValue();
		} else {
			return null;
		}
	}
	
	public String[] getProteinAlternativeNames() {
		if(entry.getProteinDescription().hasAlternativeNames()) {
			List<Name> altNames = entry.getProteinDescription().getAlternativeNames();
			String[] ans = new String[altNames.size()];
			for(int i = 0; i < altNames.size(); i++) {
				ans[i] = altNames.get(i).getFields().get(0).getValue();
			}
			return ans;
		} else {
			return null;
		}
	}
	
	public String[] getProteinSubNames() {
		if(entry.getProteinDescription().hasSubNames()) {
			List<Name> subNames = entry.getProteinDescription().getSubNames();
			String[] ans = new String[subNames.size()];
			for(int i = 0; i < subNames.size(); i++) {
				ans[i] = subNames.get(i).getFields().get(0).getValue();
			}
			return ans;
		} else {
			return null;
		}
	}
	
	public String[] getKeywords() {
		if(entry.getKeywords() != null) {
			List<Keyword> keywords = entry.getKeywords();
			String[] ans = new String[keywords.size()];
			for(int i = 0; i < keywords.size(); i++) {
				ans[i] = keywords.get(i).getValue();
			}
			return ans;
		} else {
			return null;
		}
	}
	
	public String[] getFunctionComments() {
		if(entry.getComments(CommentType.FUNCTION).size() > 0) {
			List<Comment> funComments = entry.getComments(CommentType.FUNCTION);
			String[] ans = new String[funComments.size()];
			for(int i = 0; i < funComments.size(); i++) {
				ans[i] = funComments.get(i).toString();
			}
			return ans;
		} else {
			return null;
		}
	}
	
	public String getUniprotId() {
		return entry.getUniProtId().getValue();
	}
	
	public String[] getSequenceSimilarities() {
		if(entry.getComments(CommentType.SIMILARITY).size() > 0) {
			List<Comment> seqSimComments = entry.getComments(CommentType.SIMILARITY);
			String[] ans = new String[seqSimComments.size()];
			for(int i = 0; i < seqSimComments.size(); i++) {
				ans[i] = seqSimComments.get(i).toString();
			}
			return ans;
		} else {
			return null;
		}
	}
	
	/*
	
	/**
	 * @param args
	 */
					
					/*
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		UniProtQueryUtil uq = new UniProtQueryUtil();
		BufferedReader reader = new BufferedReader(new FileReader(args[0]));
		String line;
		System.out.println("PROTEIN_ACCESSION\tUNIPROT_ID\tREC_NAME\tALT_NAME\tSUB_NAME\tORGANISM\tNCBI_TAXID\tKEYWORDS\tFUNCTIONS\tSEQ_SIMILARITY");
		while((line = reader.readLine()) != null) {
			if(uq.queryByAccessionNumber(line)) {
				System.out.println(line + "\t" + uq.getUniprotId() + "\t" + uq.getProteinRecommendedName() +
						"\t" + Arrays.toString(uq.getProteinAlternativeNames()) +
						"\t" + Arrays.toString(uq.getProteinSubNames()) +
						"\t" + uq.getOrganismScientificName() + 
						"\t" + uq.getNCBITaxId() +
						"\t" + Arrays.toString(uq.getKeywords()) + 
						"\t" + Arrays.toString(uq.getFunctionComments()) +
						"\t" + Arrays.toString(uq.getSequenceSimilarities()));
			}
		}
		reader.close();
	}

	*/

}
