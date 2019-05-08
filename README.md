# `UNM_BIOCOMP_SASA`: SaSA tools, compute molecular properties

* HDO, HAC, Rings, Molecular Complexity, Rigid bonds, polar/non-polar atom counts, ionizable atom counts, TPSA, Van der Waals volume, ABE, logP, etc.

## Compilation

* Java 1.8 and Maven 3.5 or higher are required
* Reconfigured to use ChemAxon Maven repo.
* Now compiles with ChemAxon JChem 14.7.7.0.

Compile with command:
```
$ mvn clean package
```

## Usage

Run:
```
$ java -jar unm_biocomp_sasa-0.0.1-SNAPSHOT.jar com.sunset.SaSARunner -i input.sdf -o output.sdf
```
