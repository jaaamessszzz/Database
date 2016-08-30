from Utilities.Plasmid_Utilities import Plasmid_Utilities
from db.interface import DatabaseInterface
from db.model import Feature, Plasmid_File

def add_things(asdf, GG_input):
    table_info = asdf.golden_gate_assembly(GG_input) # Part plasmids for assembly of N-nano MBP cassette from 2016713
    asdf.plasmid_checks(table_info, GG_input)
    asdf.generate_ape_file(GG_input, table_info)
    # asdf.add_cassette_plasmid_to_db(GG_input, table_info)
    asdf.add_part_plasmid_to_db(GG_input, table_info)

# Generates Gen9 WT Scaffold part plasmid input for database interactions
def generate_part_plasmids(dbi, tsession):
    Gen9_Features = tsession.query(Feature).filter(Feature.description.like("%Gen9 WT%"))

    feature_dict = {}
    for Gen9_feature in Gen9_Features:
        feature_dict[Gen9_feature.Feature_name] = {}
        feature_dict[Gen9_feature.Feature_name]['Sequence'] = Gen9_feature.Feature_sequence
        feature_dict[Gen9_feature.Feature_name]['Description'] = Gen9_feature.description

    part_types = {'WT 1F10 Chain 1/B' : ['3a', '3b'],
                  'WT 1F10 Chain 2/H' : ['3a', '3b'],
                  'WT 1F11 Chain 1/A' : ['3b'],
                  'WT 1F11 Chain 2/B' : ['3a', '3b'],
                  'WT 1F12 Chain 1/B' : ['3a'],
                  'WT 1F12 Chain 2/F' : ['3a'],
                  'WT 1G1 Chain 1/B' : ['3a', '3b'],
                  'WT 1G1 Chain 2/C' : ['3a'],
                  'WT 1G2 Chain 1/D' : ['3b'],
                  'WT 1G2 Chain 2/H' : ['3b']
                  }

    enum = 3

    for part_assembly in feature_dict:
        for part_type in part_types[part_assembly]:

            if part_type == '3a':
                GG_input = user_input('JL',
                                      ['gcatCGTCTCaAGCAGGTCTCaTATG' + feature_dict[part_assembly]['Sequence'] + 'GGTAGCGGCAGCGGCAGTTCTtGAGACCtGAGACGgcat'],
                                      'Part',
                                      'pJL' + ('0000' + str(enum))[-4:],
                                      'TEST box #1',
                                      part_type,
                                      feature_dict[part_assembly]['Description']
                                      )
                add_things(asdf, GG_input)
                enum += 1


            elif part_type == '3b':
                GG_input = user_input('JL',
                                      ['gcatCGTCTCaAGCAGGTCTCATTCT' + feature_dict[part_assembly]['Sequence'] + 'taaATCCtGAGACCtGAGACGgcat'],
                                      'Part',
                                      'pJL' + ('0000' + str(enum))[-4:],
                                      'TEST box #1',
                                      part_type,
                                      feature_dict[part_assembly]['Description']
                                      )
                add_things(asdf, GG_input)
                enum += 1
            else:
                print part_type

    database_files = tsession.query(Plasmid_File)

    for database_file in database_files:
        with open(database_file.file_name, 'w+b') as file:
            file.write(database_file.File)

def main():
    class user_input(object):
        def __init__(self, creator, input_sequences, assembly_type, UID, location, part_number=None,
                     part_description=None):
            self.creator = creator
            self.input_sequences = input_sequences
            self.assembly_type = assembly_type
            self.UID = UID
            self.location = location
            self.part_number = part_number
            self.part_description = part_description

    class make_mutations(object):
        def __init__(self, Plasmid_Feature_ID, mutation_list):
            self.Plasmid_Feature_ID = Plasmid_Feature_ID
            self.mutation_list = mutation_list

    asdf = Plasmid_Utilities()
    # asdf.primer_match('TATGCCTAGGGCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCAATGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGACTAGTGCTTGGATTCTCACCAATAAAAAACGCCCGGCGGCAACCGAGCGTTCTGAACAAATCCAGATGGAGTTCTGAGGTCATTACTGGATCTATCAACAGGAGTCCAAGCGAGCTCGATATCAAATTACGCCCCGCCCTGCCACTCATCGCAGTACTGTTGTAATTCATTAAGCATTCTGCCGACATGGAAGCCATCACAGACGGCATGATGAACCTGAATCGCCAGCGGCATCAGCACCTTGTCGCCTTGCGTATAATATTTGCCCATGGTGAAAACGGGGGCGAAGAAGTTGTCCATATTGGCCACGTTTAAATCAAAACTGGTGAAACTCACCCAGGGATTGGCTGAGACGAAAAACATATTCTCAATAAACCCTTTAGGGAAATAGGCCAGGTTTTCACCGTAACACGCCACATCTTGCGAATATATGTGTAGAAACTGCCGGAAATCGTCGTGGTATTCACTCCAGAGCGATGAAAACGTTTCAGTTTGCTCATGGAAAACGGTGTAACAAGGGTGAACACTATCCCATATCACCAGCTCACCGTCTTTCATTGCCATACGAAATTCCGGATGAGCATTCATCAGGCGGGCAAGAATGTGAATAAAGGCCGGATAAAACTTGTGCTTATTTTTCTTTACGGTCTTTAAAAAGGCCGTAATATCCAGCTGAACGGTCTGGTTATAGGTACATTGAGCAACTGACTGAAATGCCTCAAAATGTTCTTTACGATGCCATTGGGATATATCAACGGTGGTATATCCAGTGATTTTTTTCTCCATTTTAGCTTCCTTAGCTCCCTGCCCATCTCGATAACTCAAAAAATACGCCCTGTAGTGATCTTATTTCATTATGGTGAAAGTTGGAACCTCTTACGTGCCGATCAACAACGAGCAGGTCTCATTCTGAGGCCGAGCGCGTCCGAGTCTTCCACAAACAGGCCTTCGAGTACATCTCCATTGCCCTGCGCATCGATGAGGATGAGAAAGCAGGACAGAAGGAGCAAGCTGTGGAATGGTATAAGAAAGGTATTGAAGAACTGGAAAAAGGAATAGCTGTTATAGTTACAGGACAAGGTGAACAGTGTGAAAGAGCTAGACGCCTTCAAGCTAAAATGATGACTAATTTGGTTATGGCCAAGGACCGCTTACAACTTCTAGAAtaaATCCtGAGACC')

    # Create up the database session
    dbi = DatabaseInterface()
    tsession = dbi.get_session()

    GG_input = user_input('JL',
                          ['pAAG28',
                           'ASDF',
                           'pAAG54',
                           'pAAG70',
                           'pAAG22',
                           'pAAG31',
                           'pAAG24',
                           'pAAG26'
                           ],
                          'Cassette',
                          'pJL0000',
                          'BUFU'
                          )

    user_mutations = make_mutations(17, [('D', 20, 'V'),('A', 46, 'K'),('A', 141, 'F')])

    asdf.generate_mutant_sequence(tsession, user_mutations)

if __name__ == '__main__':
    main()