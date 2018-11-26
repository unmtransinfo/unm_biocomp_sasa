package com.sunset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;

public class Rows2Cols2 {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		HashSet<String> idlet = new HashSet<String>();
		HashSet<String> idcode = new HashSet<String>();
		HashMap<String, HashSet<String>> drugs1 = new HashMap<String, HashSet<String>>();
		HashMap<String, HashSet<String>> drugs2 = new HashMap<String, HashSet<String>>();
		BufferedReader reader = new BufferedReader(new FileReader(args[0]));
		String line;
		reader.readLine();
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			String drug = tokens[1];
			String tar1 = tokens[2];
			String tar2 = tokens[3];
			if(!drugs1.containsKey(drug)) {
				HashSet<String> targets1 = new HashSet<String>();
				if(tar1.length() > 3) {
					targets1.add(tar1);
					drugs1.put(drug, targets1);
				}
			} else {
				if(tar1.length() > 3) {
					drugs1.get(drug).add(tar1);
				}
			}
			if(!drugs2.containsKey(drug)) {
				HashSet<String> targets2 = new HashSet<String>();
				if(tar2.length() > 3) {
					targets2.add(tar2);
					drugs2.put(drug, targets2);
				}
			} else {
				if(tar2.length() > 3) {
					drugs2.get(drug).add(tar2);
				}
			}
			if(tar1.length() > 3) {
				idlet.add(tar1);
			}
			if(tar2.length() > 3) {
				idcode.add(tar2);
			}
		}
		reader.close();
		
		PrintStream writer = new PrintStream(args[1]);
		writer.print("SMILES\tDrugName");
		for(String tar1 : idlet) {
			writer.print("\t" + tar1);
		}
		for(String tar2 : idcode) {
			writer.print("\t" + tar2);
		}
		writer.println();
		reader = new BufferedReader(new FileReader(args[0]));
		reader.readLine();
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			HashSet<String> idlet1 = drugs1.get(tokens[1]);
			HashSet<String> idcode1 = drugs2.get(tokens[1]);
			StringBuilder sb = new StringBuilder(tokens[1]);
			for(String tar1 : idlet) {
				if(idlet1.contains(tar1)) {
					sb.append("\t1");
				} else {
					sb.append("\t0");
				}
			}
			for(String tar2 : idcode) {
				if(idcode1.contains(tar2)) {
					sb.append("\t1");
				} else {
					sb.append("\t0");
				}
			}
			writer.println(line + sb.toString());
		}
		writer.close();
		reader.close();
	}
		
}
