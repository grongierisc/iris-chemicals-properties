import pandas as pd

from sqlalchemy import create_engine,types

engine = create_engine('iris://SuperUser:SYS@localhost:62348/IRISAPP')

# read csv
df = pd.read_csv("pKaInWater.csv")

# split data into two dataframes
df_basic = df[df['basicOrAcidic'] == 'basic']
df_acidic = df[df['basicOrAcidic'] == 'acidic']

# insert into table
df_basic.to_sql('pka_basic', engine, schema="Data" ,if_exists='replace', index=False, dtype={
    'Smiles': types.VARCHAR(200),
    'pKa': types.FLOAT(8,2),
    'basicOrAcidic': types.VARCHAR(10)})

df_acidic.to_sql('pka_acidic', engine, schema="Data" ,if_exists='replace', index=False, dtype={
    'Smiles': types.VARCHAR(200),
    'pKa': types.FLOAT(8,2),
    'basicOrAcidic': types.VARCHAR(10)})
    