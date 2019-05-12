# `UNM_BIOCOMP_SASA`: SaSA tools, compute molecular properties

* HDO, HAC, Rings, Molecular Complexity, Rigid bonds, polar/non-polar atom counts, ionizable atom counts, TPSA, Van der Waals volume, ABE, logP, etc.
* Developed by Oleg Ursu

## Dependencies

* Java 1.8 and Maven 3.5 or higher are required
* ChemAxon JChem (14.7.7.0 ok).
* ChemAxon Maven Repository (requires registration and credentials).
* Oracle Maven Repository (requires registration and credentials).

## Compilation

```
mvn clean
mvn package
```

## Usage

Run:
```
$ java -jar unm_biocomp_sasa-0.0.1-SNAPSHOT.jar com.sunset.SaSARunner -i input.sdf -o output.sdf
```

## Maven configuration:

pom.xml should include:

```
    <repository>
      <id>maven.oracle.com</id>
      <releases>
        <enabled>true</enabled>
      </releases>
      <snapshots>
        <enabled>false</enabled>
      </snapshots>
      <url>https://maven.oracle.com</url>
      <layout>default</layout>
    </repository>

  </repositories>

  <pluginRepositories>
    <pluginRepository>
      <id>maven.oracle.com</id>
      <url>https://maven.oracle.com</url>
    </pluginRepository>
  </pluginRepositories>
```

settings.xml should include:

```
   <server>
      <id>maven.oracle.com</id>
      <username>ORACLE-OTN-USERNAME</username>
      <password>ORACLE-OTN-PASSWORD</password>
      <configuration>
        <basicAuthScope>
          <host>ANY</host>
          <port>ANY</port>
          <realm>OAM 11g</realm>
        </basicAuthScope>
        <httpConfiguration>
          <all>
            <params>
              <property>
                <name>http.protocol.allow-circular-redirects</name>
                <value>%b,true</value>
              </property>
            </params>
          </all>
        </httpConfiguration>
      </configuration>
    </server>
```
