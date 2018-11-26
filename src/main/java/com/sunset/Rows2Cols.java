package com.sunset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class Rows2Cols {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		HashSet<String> targets = new HashSet<String>();
		HashMap<String, HashSet<String>> drugs = new HashMap<String, HashSet<String>>();
		BufferedReader reader = new BufferedReader(new FileReader(args[0]));
		reader.readLine();
		String line;
		while((line = reader.readLine()) != null) {
			String[] tokens = line.split("\t");
			if(tokens.length < 5) {
				continue;
			}
			String drug = tokens[1];
			String target = tokens[4];
			if(!drugs.containsKey(drug)) {
				HashSet<String> drug_targets = new HashSet<String>();
				drug_targets.add(target);
				drugs.put(drug, drug_targets);
			} else {
				HashSet<String> drug_targets = drugs.get(drug);
				drug_targets.add(target);
			}
			targets.add(target);
		}
		reader.close();
		PrintStream writer = new PrintStream(args[1]);
		writer.print("molname");
		for(String target : targets) {
			writer.print("\t" + target);
		}
		writer.println();
		for(Map.Entry<String, HashSet<String>> entry : drugs.entrySet()) {
			StringBuilder sb = new StringBuilder(entry.getKey());
			for(String target : targets) {
				if(entry.getValue().contains(target)) {
					sb.append("\t1");
				} else {
					sb.append("\t0");
				}
			}
			writer.println(sb.toString());
		}
		writer.close();
	}

}
