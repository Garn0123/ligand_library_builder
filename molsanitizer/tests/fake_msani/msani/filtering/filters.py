class Filters:
    @staticmethod
    def remove_exotic_chem_to_db2(df):
        return df[~df["smiles"].str.contains("Se")]
