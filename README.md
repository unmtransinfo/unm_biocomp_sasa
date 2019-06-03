# `UNM_BIOCOMP_SASA`

Tools to compute molecular properties.

Developed by Oleg Ursu.

* HDO, HAC, Rings, Molecular Complexity, Rigid bonds, polar/non-polar atom counts, ionizable atom counts, TPSA, Van der Waals volume, ABE, logP, etc.

## Dependencies

* Java 1.8
* Maven 3.5+
* ChemAxon JChem (14.7.7.0 ok, formerly 6.2.3).
* ChemAxon Maven Repository (requires registration and credentials).
* Oracle Maven Repository (requires registration and credentials).

## Maven configuration:

* See `pom.xml` and example [`settings.xml`](doc/settings.xml).

## Compilation

```
mvn clean package
```

## Usage

Example:

```
mvn exec:java -Dexec.mainClass="com.sunset.SaSARunner" -Dexec.args="-i input.sdf -o output.sdf"
```

