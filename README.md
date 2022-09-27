# `UNM_BIOCOMP_SASA`

Tools to compute molecular properties.

Originally developed by Oleg Ursu.

* HDO, HAC, Rings, Molecular Complexity, Rigid bonds, polar/non-polar atom counts, ionizable atom counts, TPSA, Van der Waals volume, ABE, logP, SMCM, etc.

## References

* [Rapid Evaluation of Synthetic and Molecular Complexity for in Silico Chemistry, Allu and Oprea, 2005](https://pubs.acs.org/doi/10.1021/ci0501387)

## Dependencies

* Java 1.8
* Maven 3.5+
* ChemAxon JChem (14.7.7.0 ok).
* ChemAxon Maven Repository (requires registration and credentials).

## Maven configuration:

* See `pom.xml` and example [`settings.xml`](doc/settings.xml).

## Compilation

```
mvn clean package
```

## Usage

Example:

```
java -jar unm_biocomp_sasa-0.0.1-SNAPSHOT-jar-with-dependencies.jar -i input.smi -o output.tsv
```

