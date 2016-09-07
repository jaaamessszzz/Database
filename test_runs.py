from Utilities.Plasmid_Utilities import Plasmid_Utilities
from db.interface import DatabaseInterface
from db.model import Feature, Plasmid_File, Plasmid

import sys
import pprint

def add_things(asdf, input_dict, assembly_type, part_type):
    print '\n\n'
    print input_dict['sequence']
    print '\n\n'
    asdf.plasmid_checks( input_dict, assembly_type, part_type )
    asdf.add_part_plasmid_to_db(input_dict, part_type, auto_commit=False)
    
# Generates Gen9 WT Scaffold part plasmid input for database interactions
def generate_part_plasmids(asdf, dbi, tsession):
    Gen9_Features = tsession.query(Feature).filter(Feature.description.like("%Gen9 WT%"))

    feature_dict = {}
    for Gen9_feature in Gen9_Features:
        feature_dict[Gen9_feature.Feature_name] = {}
        feature_dict[Gen9_feature.Feature_name]['Sequence'] = Gen9_feature.Feature_sequence
        feature_dict[Gen9_feature.Feature_name]['Description'] = Gen9_feature.description

    import collections
    part_types_temp = (('WT 1F10 Chain 1/B' , ['3a', '3b']),
                       ('WT 1F10 Chain 2/H' , ['3a', '3b']),
                       ('WT 1F11 Chain 1/A' , ['3b']),
                       ('WT 1F11 Chain 2/B' , ['3a', '3b']),
                       ('WT 1F12 Chain 1/B' , ['3a']),
                       ('WT 1F12 Chain 2/F' , ['3a']),
                       ('WT 1G1 Chain 1/B' , ['3a', '3b']),
                       ('WT 1G1 Chain 2/C' , ['3a']),
                       ('WT 1G2 Chain 1/D' , ['3b']),
                       ('WT 1G2 Chain 2/H' , ['3b'])
                       )

    part_list = collections.OrderedDict(part_types_temp)
    enum = 17

    assembly_type = 'part'

    for part_types in part_list:
        for part_type in part_list[part_types]:
            if part_type == '3a':
                # GG_input = user_input('JL',
                #                       ['gcatCGTCTCaAGCAGGTCTCaTATG' + feature_dict[part_types]['Sequence'] + 'GGTAGCGGCAGCGGCAGTTCTtGAGACCtGAGACGgcat'],
                #                       'Part',
                #                       'pJL' + ('0000' + str(enum))[-4:],
                #                       'TEST box #1',
                #                       part_type,
                #                       feature_dict[part_types]['Description']
                #                       )

                input_sequences = ['gcatCGTCTCaAGCAGGTCTCaTATG' + feature_dict[part_types]['Sequence'] + 'GGTAGCGGCAGCGGCAGTTCTtGAGACCtGAGACGgcat']
                table_info = asdf.golden_gate_assembly(input_sequences, assembly_type)

                input_dict = {'creator': 'JL',
                              'plasmid_name': 'pJL' + ('0000' + str(enum))[-4:],
                              'plasmid_type': assembly_type,
                              'location': 'Mah box',
                              'description': table_info[
                                                 'Complete Description'] + 'Part plasmids for assembly of N-nano MBP cassette from 2016713',
                              'sequence': table_info['Complete Assembly'],
                              'status': 'designed'
                              }

                add_things(asdf, input_dict, assembly_type, part_type)
                enum += 1


            elif part_type == '3b':
                # GG_input = user_input('JL',
                #                       ['gcatCGTCTCaAGCAGGTCTCATTCT' + feature_dict[part_types]['Sequence'] + 'taaATCCtGAGACCtGAGACGgcat'],
                #                       'Part',
                #                       'pJL' + ('0000' + str(enum))[-4:],
                #                       'TEST box #1',
                #                       part_type,
                #                       feature_dict[part_types]['Description']
                #                       )
                input_sequences = ['gcatCGTCTCaAGCAGGTCTCATTCT' + feature_dict[part_types]['Sequence'] + 'taaATCCtGAGACCtGAGACGgcat']
                table_info = asdf.golden_gate_assembly(input_sequences, assembly_type)

                input_dict = {'creator': 'JL',
                              'plasmid_name': 'pJL' + ('0000' + str(enum))[-4:],
                              'plasmid_type': assembly_type,
                              'location': 'Mah box',
                              'description': table_info[
                                                 'Complete Description'] + 'Part plasmids for assembly of N-nano MBP cassette from 2016713',
                              'sequence': table_info['Complete Assembly'],
                              'status': 'designed'
                              }

                add_things(asdf, input_dict, assembly_type, part_type)
                enum += 1
            else:
                print part_type

    database_files = tsession.query(Plasmid_File)

    for database_file in database_files:
        with open(database_file.file_name, 'w+b') as file:
            file.write(database_file.File)

def generate_cassette_plasmids(asdf, dbi, tsession):
    input_sequences = ['pAAG28',
                       'pAAG67',
                       'pJL0003',
                       'pAAG52',
                       'pAAG22',
                       'pAAG31',
                       'pAAG87'
                       ]

    assembly_type = 'cassette'

    table_info = asdf.golden_gate_assembly(input_sequences, assembly_type)

    input_dict = {'creator': 'JL',
                  'plasmid_name': 'pJL0017',
                  'plasmid_type': assembly_type,
                  'location': 'Mah box',
                  'description': table_info[
                                     'Complete Description'],
                  'sequence': table_info['Complete Assembly'],
                  'status': 'designed'
                  }

    asdf.plasmid_checks(input_dict, assembly_type)
    asdf.add_cassette_plasmid_to_db(input_dict, table_info, auto_commit=True)

def main():

    class make_mutations(object):
        def __init__(self, Plasmid_Feature_ID, mutation_list):
            self.Plasmid_Feature_ID = Plasmid_Feature_ID
            self.mutation_list = mutation_list

    asdf = Plasmid_Utilities()

    # Create up the database session
    dbi = DatabaseInterface()
    tsession = dbi.get_session()


    asdf.generate_ape_from_database_ID(tsession, 'AG', '3', write_to_file=True)
    # generate_cassette_plasmids(asdf, dbi, tsession)
    # generate_part_plasmids(asdf, dbi, tsession)

    # generate_cassette_plasmids(asdf, dbi, tsession)
    # generate_part_plasmids(asdf, dbi, tsession)

    # user_mutations = make_mutations(306, [('D', 20, 'V'),('A', 46, 'K'),('A', 141, 'F')]) # 1G2 Chain 1
    # asdf.generate_mutant_sequence(tsession, user_mutations)

if __name__ == '__main__':
    main()
