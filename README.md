# SaSA tools, compute molecular properties

- HDO, HAC, Rings, Molecular Complexity, Rigid bonds, polar/non-polar atom counts, ionizable atom counts, TPSA, Van der Waals volume, ABE, logP, etc.

## Compilation
- Java 1.8 and Maven 3.5 or higher are required
- * Edit pom.xml file and rename iridium repo to localhost and change location to point to unzipped maven2.zip in /home/data/oleg/maven2.zip

Compile with command:
```
$ mvn clean package
```

## Usage
Run:
```
$ java -jar SaSA-1.0.jar com.sunset.SaSARunner -i input.sdf -o output.sdf
```