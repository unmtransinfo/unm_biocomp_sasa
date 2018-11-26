package com.sunset;

import java.util.HashMap;
import java.util.Map;

public enum VDW {
	HYDROGEN(1, "H", 1.2d, 7.24d),
	LITHIUM(2, "Li", 1.82d, 25.25d),
	HELIUM(2, "He", 1.4d, 11.49d),
	SODIUM(11, "Na", 2.27d, 49.00d),
	MAGNESIUM(12, "Mg", 1.73d, 21.69d),
	POTASSIUM(19, "K", 2.75d, 87.11d),
	NICKEL(28, "Ni", 1.63d, 18.14d),
	COPPER(29, "Cu", 1.4d, 11.49d),
	ZINC(30, "Zn", 1.39d, 11.25d),
	GALLIUM(31, "Ga", 1.87d, 27.39d),
	PALLADIUM(46, "Pd", 1.63d, 18.14d),
	SILVER(47, "Ag", 1.72d, 21.31d),
	CADMIUIM(48, "Cd", 1.58d, 16.52d),
	INDIUM(49, "In", 1.93d, 30.11d),
	TIN(50, "Sn", 2.17d, 42.80d),
	PLATINUM(78, "Pt", 1.72d, 21.31d),
	GOLD(79, "Au", 1.66d, 19.16d),
	THALLIUM(81, "Tl", 1.96d, 31.54d),
	LEAD(82, "Pb", 2.02d, 34.53d),
	URANIUM(92, "U", 1.86d, 26.95d),
	BERYLLIUM(4, "Be", 1.53d, 15.00d),
	BORON(5, "B", 1.92d, 29.65d),
	ALUMINIUM(13, "Al", 1.84d, 26.09d),
	CALCIUM(20, "Ca", 2.31d, 51.63d),
	GERMANIUM(32, "Ge", 2.11d, 39.35d),
	RUBIDIUM(37, "Rb", 3.03d, 116.52d),
	STRONTIUM(38, "Sr", 2.49d, 64.67d),
	ANTIMONY(51, "Sb", 2.06d, 36.62d),
	CAESIUM(55, "Cs", 3.43d, 169.03d),
	BARIUM(56, "Ba", 2.68d, 80.63d),
	BISMUTH(83, "Bi", 2.07d, 37.15d),
	POLONIUM(84, "Po", 1.97d, 32.02d),
	ASTATINE(85, "At", 2.02d, 34.53d),
	RADON(86, "Rn", 2.2d, 44.60d),
	FRANCIUM(87, "Fr", 3.48d, 176.53d),
	RADIUM(88, "Ra", 2.83d, 94.94d),
	CARBON(6, "C", 1.7d, 20.58),
	NITROGEN(7, "N", 1.55d, 15.6d),
	OXYGEN(8, "O", 1.52d, 14.71d),
	FLUORINE(9, "F", 1.47d, 13.31d),
	NEON(10, "Ne", 1.54d, 15.30d),
	CHLORINE(17, "Cl", 1.75d, 22.45d),
	ARGON(18, "Ar", 1.88d, 27.83d),
	BROMINE(35, "Br", 1.85d, 26.52d),
	IODINE(53, "I", 1.98d, 32.52d),
	PHOSPHORUS(15, "P", 1.80d, 24.43d),
	SULFUR(16, "S", 1.80d, 24.43d),
	ASRSENIC(33, "As", 1.85d, 26.52d),
	SILICON(14, "Si", 2.10d, 38.79d),
	SELENIUM(34, "Se", 1.90d, 28.73d),
	TELLURIUM(52, "Te", 2.06d, 36.62d),
	XENON(54, "Xe", 2.16d, 42.21d);

	private final int number;
	
	private final String symbol;
	
	private final double radius;
	
	private final double volume;
	
	static final Map<String, VDW> symbolMap = new HashMap<>();
	
	static final Map<Integer, VDW> numberMap = new HashMap<>();
	
	static {
		for(final VDW value : values()) {
			numberMap.put(value.number, value);
			symbolMap.put(value.symbol, value);
		}
	}
	
	private VDW(int number, String symbol, double radius, double volume) {
		this.number = number;
		this.symbol = symbol;
		this.radius = radius;
		this.volume = volume;
	}
	
	public int number() {
		return number;
	}
	
	public String symbol() {
		return symbol;
	}
	
	public double radius() {
		return radius;
	}
	
	public double volume() {
		return volume;
	}
	
	public static VDW ofNumber(final int number) {
		return numberMap.get(number);
	}
	
	public static VDW ofSymbol(final String symbol) {
		return symbolMap.get(symbol);
	}
	
}
