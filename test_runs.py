from Utilities.Plasmid_Utilities import Plasmid_Utilities
from Utilities.Plasmid_View_Tools import Plasmid_View_Tools
from db.interface import DatabaseInterface
from db.model import Feature, Plasmid_File, Plasmid, Part_Plasmid_Part, CDS_Mutant_Constituent
from sqlalchemy import and_

import sys
import pprint
import traceback

def add_things(asdf, input_dict, assembly_type, part_type):
    print '\n\n'
    print input_dict['sequence']
    print '\n\n'
    asdf.plasmid_checks( input_dict, assembly_type, part_type )
    asdf.add_part_plasmid_to_db(input_dict, part_type, auto_commit=False)
    
# Generates Gen9 WT Scaffold part plasmid input for database interactions
def generate_part_plasmids(asdf, tsession):
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
    enum = 17 # Start UIDs at index 17

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

def generate_cassette_plasmids(asdf):

    input_sequences = [('AG', 1),
                       ('AG', 3),
                       ('AG', 15),
                       ('AG', 19),
                       ('AG', 8),
                       ('AG', 10),
                       ('AG', 17)
                       ]

    assembly_type = 'cassette'

    table_info = asdf.golden_gate_assembly(input_sequences, assembly_type)

    input_dict = {'creator': 'JL',
                  'plasmid_name': 'pJLTEST',
                  'plasmid_type': assembly_type,
                  'location': 'Mah box',
                  'description': table_info[
                                     'Complete Description'],
                  'sequence': table_info['Complete Assembly'],
                  'status': 'designed'
                  }

    asdf.plasmid_checks(input_dict, assembly_type)

    import pprint
    pprint.pprint(input_dict['sequence'])

    # asdf.add_cassette_plasmid_to_db(input_dict, table_info, auto_commit=False)

def generate_multicassette_plasmids(asdf):

    assembly_type = 'multicassette'
    input_sequences = [('JL', 28), ('JL', 30), ('JL', 44)]

    table_info = asdf.golden_gate_assembly(input_sequences, assembly_type)

    input_dict = {'creator': 'JL',
                  'plasmid_name': 'pJLTEST2',
                  'plasmid_type': assembly_type,
                  'location': 'Mah box',
                  'description': table_info['Complete Description'],
                  'sequence': table_info['Complete Assembly'],
                  'status': 'designed'
                  }

    asdf.plasmid_checks(input_dict, assembly_type)

    import pprint
    pprint.pprint(input_dict['sequence'])

    asdf.add_multicassette_plasmid_to_db(input_dict, table_info, auto_commit=True)

def main():

    # Create up the database session
    dbi = DatabaseInterface()
    try:
        tsession = dbi.get_session()

        asdf = Plasmid_Utilities(tsession)
        qwer = Plasmid_View_Tools(tsession)

        # asdf.generate_ape_from_database_ID('JL', '1', write_to_file=True)
        # generate_part_plasmids(asdf)
        # generate_cassette_plasmids(asdf)
        # generate_multicassette_plasmids(asdf)
        #asdf.update_features()
        # asdf.nuke_plasmids([('JL', 59), ('JL', 60), ('JL', 61), ('JL', 62), ('JL', 63)])

        # feature_indicies_list, part_indicies_list = qwer.get_plasmid_indicies(('AG', 3))
        # print feature_indicies_list
        # print part_indicies_list

        ### Generate Primers for a given sequence and TM
        # target_primer_F, target_primer_R = qwer.generate_primers('ggcgtcatgactgatgtccaccggcgcttcctccagttgctgatgacgcatggcgtgctagaggaatgggacgtgaagcgcttgcagacgcactgctacaaggtccatgaccgcaatgccaccgtagataagttggaggacttcatcaacaacattaacagtgtcttggagtccttgtatattgagataaagagaggagtcacggaagatgatggcagacccatttatgcgttggtgaatcttgctacaacttcaatttccaaaatggctacggattttgcagagaatgaactggatttgtttagaaaggctctggaactgattattgactcagaaaccggctttgcgtcttccacaaacatattgaacctggttgatcaacttaaaggcaagaagatgaggaagaaggaagcggagcaggtgctgcagaagtttgttcaaaacaagtggctgattgagaaggaaggggagttcaccctgcacggccgggccatcctggagatggagcaatacatccgggaaacgtaccccgacgcggtgaagatctgcaatatctgtcacagcctcctcatccagggtcaaagctgcgaaacctgtgggatcaggatgcacttaccctgcgtggccaagtacttccagtcgaatgctgaaccgcgctgcccccactgcaacgactactggccccacgagatcccaaaagtcttcgaccctgag',
        #                                                          60,
        #                                                          left_arm='gcatCGTCTCaAGCAGGTCTCATTCT', # 5' -> 3'
        #                                                          right_arm='taaATCCtGAGACCtGAGACGgcat'  # 5' -> 3'
        #                                                          )
        # With defaults at TM = 60
        # F Generated: GCATCGTCTCAAGCAGGTCTCATTCTGGCGTCATGACTGATGTCCACC
        # F Designed : gcatCGTCTCaAGCAGGTCTCATTCTGGCGTCATGACTGATGTCCACC
        # R Generated: ATGCCGTCTCAGGTCTCAGGATTTACTCAGGGTCGAAGACTTTTGGGA
        # R Designed : atgcCGTCTCaGGTCTCaGGATttaCTCAGGGTCGAAGACTTTTGGGAT

        ### Generate a Mutant Plasmid Feature and Associated .ape file, push to database
        # Plasmid_Feature_ID = 984 # WT 1F10 Chain 1/B
        # mutation_list = [('L', 20, 'R'),('K', 40, 'A'),('C', 60, 'D')]
        # Mutant_plasmid_sequence, CDS_mutant_constituents, mutant_genbank_file = asdf.generate_mutant_sequence(Plasmid_Feature_ID, mutation_list)
        # asdf.add_mutant_to_db( Plasmid_Feature_ID, CDS_mutant_constituents, auto_commit=False )

        ### Feature Design!
        # input_dict = [{'feature_ID': 2862, 'design_sequence': '|---T7 Promoter---|', 'design_name': 'ASDFASDF', 'design_description': 'HERP DERP'},
        #               {'feature_ID': 2860, 'design_sequence': '|---Lac Operon---|', 'design_name': 'ASDFASDF', 'design_description': 'HERP DERP'},
        #               {'feature_ID': 2856, 'design_sequence': '|---LgBiT---|', 'design_name': 'ASDFASDF', 'design_description': 'HERP DERP'},
        #               {'feature_ID': 2857, 'design_sequence': '|---SmBiT---|', 'design_name': 'ASDFASDF', 'design_description': 'HERP DERP'},
        #               {'feature_ID': 2858, 'design_sequence': '|---WT 1F10 Chain 1/B---|', 'design_name': 'ASDFASDF', 'design_description': 'HERP DERP'},
        #               {'feature_ID': 2859, 'design_sequence': '|---WT 1F10 Chain 2/H---|', 'design_name': 'ASDFASDF', 'design_description': 'HERP DERP'},
        #               {'feature_ID': 2861, 'design_sequence': '|---ColE1---|', 'design_name': 'ASDFASDF', 'design_description': 'HERP DERP'},
        #               {'feature_ID': 2863, 'design_sequence': '|---KAN---|', 'design_name': 'ASDFASDF', 'design_description': 'HERP DERP'},
        #               {'feature_ID': 2864, 'design_sequence': '|---Terminator---|', 'design_name': 'ASDFASDF', 'design_description': 'HERP DERP'}
        #               ]
        # asdfasdf = asdf.design_feature(input_dict, 'pJLASDFASDF', 'ALL YOUR DERPS BELONG TO US')

        input_dict = [{'feature_ID': 3031, 'design_sequence': 'GTTTCAGAGCATGCTGGAAACAGCATAGCAAGTTGAAATAAGGTCTTCCCCGCATCCGCCGATACCAGCCGAAAGGCCCTTGGCAGCGACGGCACCGAGTCGGTGCTTTTTT', 'design_name': 'mhf/37', 'design_description': 'Some design from Kale to test.'}]
        asdfasdf = asdf.design_feature(input_dict, 'pJL_Kale-test', 'Test design from Kale', 'ASDF')

        # my_plasmids = tsession.query(Plasmid, Part_Plasmid_Part).filter(and_(Plasmid.creator == Part_Plasmid_Part.creator, Plasmid.creator_entry_number == Part_Plasmid_Part.creator_entry_number))
        # for plasmid, part_plasmid_part in my_plasmids:
        #     input_dict = {'creator': plasmid.creator,
        #                   'plasmid_name': plasmid.plasmid_name,
        #                   'plasmid_type': plasmid.plasmid_type,
        #                   'location': plasmid.location,
        #                   'description': plasmid.description,
        #                   'sequence': plasmid.sequence,
        #                   'status': plasmid.status
        #                   }
        #     try:
        #         if plasmid.plasmid_type == 'part':
        #             asdf.plasmid_checks(input_dict, plasmid.plasmid_type, part_plasmid_part.part_number)
        #         else:
        #             asdf.plasmid_checks(input_dict, plasmid.plasmid_type)
        #     except:
        #         print '{0}{1}'.format(plasmid.creator, plasmid.creator_entry_number)

        # print(asdf.return_designable_features(dbi.get_engine(), 'TP', 1))

        tsession.close()
    except Exception, e:
        del dbi
        print(str(e))
        print(traceback.format_exc())

if __name__ == '__main__':
    main()

