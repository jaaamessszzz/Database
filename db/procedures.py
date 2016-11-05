stored_procedure_column_names = dict(
    getPartPlasmidParts = [
        'Plasmid_creator', 'Plasmid_creator_entry_number', 'Plasmid_plasmid_name', 'Plasmid_plasmid_type', 'Plasmid_location', 'Plasmid_description', 'Plasmid_sequence', 'Plasmid_sequence_UC', 'Plasmid_status', 'Plasmid_date',
        'Part_Plasmid_creator_entry_number', 'Part_Plasmid_creator', 'Part_Plasmid_resistance',
        'Part_Plasmid_Part_creator_entry_number', 'Part_Plasmid_Part_creator', 'Part_Plasmid_Part_part_number',
        'Part_Type_part_number', 'Part_Type_overhang_5', 'Part_Type_overhang_3'],
    getDesignableFeatures = ['feature_ID', 'creator', 'creator_entry_number', 'feature_name', 'MD5_hash', 'sequence', 'description', 'feature_type', 'color'],
)

def call_procedure(engine, function_name, params, as_dict = False):
    connection = engine.raw_connection()
    try:
        cursor = connection.cursor()
        cursor.callproc(function_name, params)
        results = cursor.fetchall()

        # Convert the results to a dictionary if requested. This requires the query columns to be properly defined in stored_procedure_column_names so it is fragile.
        if as_dict:
            assert(function_name in stored_procedure_column_names)
            column_names = stored_procedure_column_names[function_name]
            if results:
                dresults = []
                assert(len(results[0]) == len(column_names))
                results = [dict(zip(column_names, r)) for r in results]

        cursor.close()
        connection.commit()
        return results
    finally:
        connection.close()



