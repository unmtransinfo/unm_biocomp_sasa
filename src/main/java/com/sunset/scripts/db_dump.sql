-- export data tables from WBPK except STRUCTURES table

connect 'jdbc:derby:db';
AUTOCOMMIT OFF;
CALL SYSCS_UTIL.SYSCS_EXPORT_QUERY('SELECT STRUCT_ID,SWISSPROTID,PDBID,TARGET_NAME,TARGET_CLASS,TARGET_BIOACTIVITY_VALUE,TARGET_BIOACTIVITY_TYPE,CONFIDENCE,DATA_SOURCE,ACTVSMEDI,ACTVSTSC,ACTVSBOTH FROM DRUG_TARGET','DRUG_TARGET.tab','	',NULL,NULL);
CALL SYSCS_UTIL.SYSCS_EXPORT_QUERY('SELECT STRUCT_ID,REF,ACT_TYPE,ACT_VAL FROM HERG_BIOACT','HERG_BIOACT.tab','	',NULL,NULL);
CALL SYSCS_UTIL.SYSCS_EXPORT_QUERY('SELECT STRUCT_ID,MIXTURE_ID,FILTER,DERIVATIVE,TRADE_NAME,CASNO,DRUG_CODE,MELTING_POINT FROM MIXTURE','MIXTURE.tab','	',NULL,NULL);
CALL SYSCS_UTIL.SYSCS_EXPORT_QUERY('SELECT MIXTURE_ID,VALUE,UNIT,TYPE,SPECIES,COMMENT FROM MIXTURE_LD50','MIXTURE_LD50.tab','	',NULL,NULL);
CALL SYSCS_UTIL.SYSCS_EXPORT_QUERY('SELECT MIXTURE_ID,VALUE,UNIT,COMMENT FROM MIXTURE_WATER_SOL','MIXTURE_WATER_SOL.tab','	',NULL,NULL);
CALL SYSCS_UTIL.SYSCS_EXPORT_QUERY('SELECT STRUCT_ID,NAME,SWISSPROTID FROM PHASE1_METAB_ENZYME','PHASE1_METAB_ENZYME.tab','	',NULL,NULL);
CALL SYSCS_UTIL.SYSCS_EXPORT_QUERY('SELECT STRUCT_ID,REF,TARGET_TYPE,TARGET_NAME,ACT_TYPE,ACT_VAL,INHIB_PERC,LIG_EFF FROM WOMBAT_ACTLIST','WOMBAT_ACTLIST.tab','	',NULL,NULL);
CALL SYSCS_UTIL.SYSCS_EXPORT_QUERY('SELECT STRUCT_ID,TRADE_NAME,SYEAR,VALUE,REPORT_COMPANY FROM SALES','SALES.tab','	',NULL,NULL);
QUIT;