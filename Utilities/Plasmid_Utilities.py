import sys
import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC

sys.path.insert(0, '..')
try:
    from db.interface import DatabaseInterface
    from db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature, Plasmid_Feature, Plasmid_File
except:
    # nasty hack since we are not packaging things up properly yet for external use (e.g. the website)
    from kprimers.db.interface import DatabaseInterface
    from kprimers.db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature, Plasmid_Feature, Plasmid_File


class Plasmid_Utilities(object):
    '''
    This is a class with modules for completing tasks involving manipulation of Plasmid sequences.

    ###########################
    # Available Funtionality: #
    ###########################

    * Find available cloning primers for a plasmid of interest
    * Golden Gate Assembly within a modular cloning (MoClo) system - see http://pubs.acs.org/doi/abs/10.1021/sb500366v
    * Verify that Plasmids adhere to MoClo sequence standards
    * Generate annotated .ape files based on desired sequence features
    * Upload assembled plasmids to the database (part and cassette plasmids implemented)

    ###########################
    # Planned Funtionality:   #
    ###########################

    * Multicassette assembly and submission
    * Annotating working and archive locations for all plasmids
    * Storing information for mutations and libraries
    '''
    def __init__(self, tsession = None):
        # Create up the database session
        self.dbi = None
        try:
            self.dbi = DatabaseInterface()
        except:
            if not tsession:
                raise
        self.tsession = tsession or self.dbi.get_session()

    def reverse_complement(self, input_sequence):
        base_pair = {'A': 'T',
                     'T': 'A',
                     'C': 'G',
                     'G': 'C'
                     }
        return ''.join(reversed([base_pair[base.upper()] for base in input_sequence]))

    def primer_match(self, target_sequence):
        '''
        This function will take a target plasmid and search the Primers database for compatible cloning primers.

        :param target_sequence: sequence of the target plasmid
        :return:
        '''

        # Create up the database session
        dbi = DatabaseInterface()
        tsession = dbi.get_session()

        # Create a map from usernames to the database IDs (typically initials)
        user_map = {}
        for u in tsession.query(Users):
            user_map[u.lab_username] = u.ID

        primer_sequences = tsession.query(Primers.creator, Primers.creator_entry_number, Primers.sequence)

        target_reverse_complement = self.reverse_complement(target_sequence).upper()

        # Added [:20] since I don't currently have any cloning primers in the database
        for creator, entry_number, sequence in primer_sequences:
            ID = (creator, entry_number)

            if sequence.upper()[20:] in target_sequence.upper():
                print sequence[20:].upper()
                print 'Primer', ID, "(F) binds the target plasmid at position", target_sequence.upper().find(sequence.upper()[20:])
                print '\n'

            if sequence.upper()[20:] in target_reverse_complement:
                print sequence[20:].upper()
                print 'Primer', ID, "(R) binds the target plasmid with 5' at position", (len(target_sequence) - target_reverse_complement.find(sequence.upper()[20:]))
                print '\n'

    def golden_gate_assembly(self, user_input):
        '''
        This function will take user input of a list of plasmids and attempt to perform a golden gate assembly.

        class user_input(object):
            def __init__(self, creator, input_sequences, assembly_type, UID, location):
                self.creator = creator (creator initials)
                self.input_sequences = input_sequences (user defined names of input plasmids for assembly)
                self.assembly_type = assembly_type (type of assembly e.g. part, cassette, multicassette)
                self.UID = UID (user defined name for the assembled plasmid)
                self.location = location (where is the user planning on putting this new plasmid)

        :return: complete plasmid assembly as string
        '''

        table_info = {}
        table_info['Part list'] = []
        table_info['Part list ID'] = []
        table_info['Description'] = []
        table_info['Connectors'] = {}

        def fetch_cassette_parts():
            # Get sequences and part numbers for listed part plasmids
            part_plasmids_query = self.tsession.query(Plasmid, Part_Plasmid, Part_Plasmid_Part, Part_Type)\
                .filter(Part_Plasmid.creator == Plasmid.creator)\
                .filter(Part_Plasmid.creator_entry_number == Plasmid.creator_entry_number)\
                .filter(Part_Plasmid.creator == Part_Plasmid_Part.creator)\
                .filter(Part_Plasmid.creator_entry_number == Part_Plasmid_Part.creator_entry_number)\
                .filter(Plasmid.plasmid_name.in_(user_input.input_sequences)) \
                .filter(Part_Plasmid_Part.part_number == Part_Type.part_number)

            restriction_enzyme = 'BsaI'
            site_F = 'GGTCTC'
            site_R = 'GAGACC'

            for plasmid, part_plasmid, part_plasmid_part, part_type in part_plasmids_query:

                print "%s\t%s\t%s\t%s"% (plasmid.plasmid_name, plasmid.creator, plasmid.creator_entry_number, part_plasmid_part.part_number)

                assert plasmid.sequence.count(site_F) == 1,\
                    'There is more than one forward %s site in %s!' % (restriction_enzyme, plasmid.plasmid_name)
                assert plasmid.sequence.count(site_R) == 1, \
                    'There is more than one reverse %s site in %s!' % (restriction_enzyme, plasmid.plasmid_name)

                sequence_upper = plasmid.sequence.upper()
                site_F_position = sequence_upper.find(part_type.overhang_5, sequence_upper.find(site_F))
                site_R_position = sequence_upper.find(part_type.overhang_3, sequence_upper.find(site_R) - 10) + 4

                # Circular permutation so that forward cut site is alway upstream of reverse cut site in linear sequence
                if site_F_position > site_R_position:
                    sequence_upper = sequence_upper[site_R_position + 7:] + sequence_upper[:site_R_position + 7]

                table_info['Part list'].append( sequence_upper[
                                  sequence_upper.find(part_type.overhang_5, sequence_upper.find(site_F)):
                                  sequence_upper.find(part_type.overhang_3, sequence_upper.find(site_R) - 10) + 4
                                  ])
                table_info['Description'].append( plasmid.description )
                table_info['Part list ID'].append ((part_plasmid_part.part_number,plasmid.creator, plasmid.creator_entry_number))
                if part_plasmid_part.part_number == '1' or part_plasmid_part.part_number == '5':
                    table_info['Connectors'][part_plasmid_part.part_number] = [plasmid.creator, plasmid.creator_entry_number]

            return table_info

        # Start with Part one (cassette assembly) so that there is consistency in how the plasmids are assembled
        if user_input.assembly_type.lower() == 'cassette':
            table_info = fetch_cassette_parts()
            for part in table_info['Part list']:
                if part[:4].upper() == 'CCCT':
                    intermediate = part.upper()
                    table_info['Part list'].remove(part)
                    break

        if user_input.assembly_type.lower() == 'part':
            part_entry_vector = self.tsession.query(Plasmid).filter(Plasmid.creator == 'JL').filter(Plasmid.creator_entry_number == 2)
            sequence_upper = user_input.input_sequences[0].upper()
            table_info['Part list'].append(sequence_upper[ sequence_upper.find('CGTCTC') + 7 : sequence_upper.find('GAGACG') - 1])
            for asdf in part_entry_vector:
                intermediate = asdf.sequence[ : (asdf.sequence.upper().find('GAGACG') - 1) ]
                print intermediate

        # if user_input.assembly_type == 'Multicassette':
        #     site_F = 'CGTCTC'
        #     site_R = 'GAGACG'
        #
        #     for part in part_list:
        #         if part[:4].upper() == 'CTGA':
        #             intermediate = part.upper()
        #             table_info['Part list'].remove(part)
        #             break

        assembly_step = 0
        assembly_steps = len(table_info['Part list'])

        while assembly_step < assembly_steps:
            for part in table_info['Part list']:
                # If part 5' and intermediate 3' match...
                if part.upper()[:4] == intermediate.upper()[-4:]:
                    intermediate = intermediate.upper()[:-4] + part.upper()
            assembly_step += 1

        assert intermediate[-4:] == intermediate[:4], 'Incomplete assembly! This assembly does not produce a circular plasmid! :('

        table_info['Complete Assembly'] = intermediate[4:].upper()
        table_info['Complete Description'] = '\n'.join(table_info['Description'])

        return table_info

    def plasmid_checks(self, table_info, user_input):
        '''
        Verifies that a plasmid conforms to modular cloning sequence standards before entry into database

        # Part Plasmids #
          * Part 1 and Part 5 should contain one reverse BsmBI site
          * If not a Part 1 or Part 5, sequence should not contain BsmBI sites
          *
        :return:
        '''
        if user_input.assembly_type.lower() == 'part':
            if user_input.part_number == '1' or user_input.part_number == '5':
                assert table_info['Complete Assembly'].count('GAGACG') == 1, \
                    'There should be only one reverse %s site in a Connector Part! Your part %s submission contains %s reverse %s sites.' % (
                    'BsmBI', user_input.part_number, table_info['Complete Assembly'].count('GAGACG'), 'BsmBI')
            else:
                assert table_info['Complete Assembly'].count('GAGACG') == 0, \
                    'There should not be any reverse %s sites in this part type! Your submission %s contains %s reverse %s site(s).' % (
                        'BsmBI', user_input.UID, table_info['Complete Assembly'].count('GAGACG'), 'BsmBI')
            assert table_info['Complete Assembly'].count('CGTCTC') == 0, \
                'There is more than one forward %s site in %s! This plasmid contains %s forward %s sites.' % (
                    'BsmBI', user_input.UID, table_info['Complete Assembly'].count('CGTCTC'), 'BsmBI')
            # assert table_info['Complete Assembly'].count('GGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGACTAGTGCTTGGATTCTCACCAATAAAAAACGCCCGGCGGCAACCGAGCGTTCTGAACAAATCCAGATGGAGTTCTGAGGTCATTACTGGATCTATCAACAGGAGTCCAAGCGAGCTCGATATCAAATTACGCCCCGCCCTGCCACTCATCGCAGTACTGTTGTAATTCATTAAGCATTCTGCCGACATGGAAGCCATCACAAACGGCATGATGAACCTGAATCGCCAGCGGCATCAGCACCTTGTCGCCTTGCGTATAATATTTGCCCATGGTGAAAACGGGGGCGAA') == 0, \
            #     "You do not need to include the backbone sequence in your submission, we'll take care of that!" # 607-1108 of pAAG16, it's just a rondom section of the backbone that overlaps both the origin and resistance marker
            assert table_info['Complete Assembly'].count('GGTCTC') == 1, \
                'There is more than one forward %s site in %s! Your submission contains %s forward %s sites.' % (
                    'BsaI', user_input.UID, table_info['Complete Assembly'].count('GGTCTC'), 'BsaI')
            assert table_info['Complete Assembly'].count('GAGACC') == 1, \
                'There is more than one reverse %s site in %s! Your submission contains %s reverse %s sites.' % (
                    'BsaI', user_input.UID, table_info['Complete Assembly'].count('GAGACC'), 'BsaI')


        if user_input.assembly_type.lower() == 'cassette':
            assert table_info['Complete Assembly'].count('CGTCTC') == 1, \
                'There is more than one forward %s site in %s! Your submission contains %s forward %s sites.' % ('BsmBI', user_input.UID, table_info['Complete Assembly'].count('CGTCTC'), 'BsmBI')
            assert table_info['Complete Assembly'].count('GAGACG') == 1, \
                'There is more than one reverse %s site in %s! Your submission contains %s reverse %s sites.' % ('BsmBI', user_input.UID, table_info['Complete Assembly'].count('GAGACG'), 'BsmBI')

        if user_input.assembly_type.lower() == 'multicassette':
            # Ehhh... I'll leave this for now. There really aren't any restrictions for a multicasssette plasmid that I can think of for now
            assert table_info['Complete Assembly'].count('GGTCTC') == 1, \
                'There is more than one forward %s site in %s! Your submission contains %s forward %s sites.' % (
                'BsaI', user_input.UID, table_info['Complete Assembly'].count('GGTCTC'), 'BsaI')
            assert table_info['Complete Assembly'].count('GAGACC') == 1, \
                'There is more than one reverse %s site in %s! Your submission contains %s reverse %s sites.' % (
                'BsaI', user_input.UID, table_info['Complete Assembly'].count('GAGACC'), 'BsaI')


    def add_part_plasmid_to_db(self, user_input, table_info, auto_commit = True):
        new_plasmid_entry = {'creator': user_input.creator,
                             'plasmid_name': user_input.UID,
                             'plasmid_type': user_input.assembly_type,
                             'location': user_input.location,
                             'description': user_input.part_description,
                             'sequence': table_info['Complete Assembly'],
                             'status': 'designed'
                             }

        current_plasmid_entry = Plasmid.add(self.tsession, new_plasmid_entry, silent=True)
        user_input.database_ID = Plasmid.get_id()

        new_part_plasmid_entry = {'creator' : user_input.creator,
                                  'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                  'resistance' : 'CM'
                                  }
        Part_Plasmid.add(self.tsession, new_part_plasmid_entry, silent=True)

        for part in user_input.part_number.split(','):
            new_part_plasmid_part_entry = {'creator' : user_input.creator,
                                           'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                           'part_number' : part.strip()
                                           }
            Part_Plasmid_Part.add(self.tsession, new_part_plasmid_part_entry, silent = True)

        self.add_features(table_info, current_plasmid_entry)
        self.upload_file(user_input, table_info, current_plasmid_entry)

        if auto_commit:
            self.tsession.commit()

    def add_cassette_plasmid_to_db(self, user_input, table_info):
        new_plasmid_entry = {'creator' : user_input.creator,
                             'plasmid_name' : user_input.UID,
                             'plasmid_type' : user_input.assembly_type,
                             'location' : user_input.location,
                             'description' : table_info['Complete Description'],
                             'sequence' : table_info['Complete Assembly'],
                             'status' : 'designed'
                             }

        current_plasmid_entry = Plasmid.add(self.tsession, new_plasmid_entry, silent=False)

        from sqlalchemy import and_
        left_connector = self.tsession.query(Cassette_Connector)\
            .filter(and_(Cassette_Connector.creator == table_info['Connectors']['1'][0],
                         Cassette_Connector.creator_entry_number == table_info['Connectors']['1'][1])
                    )
        right_connector = self.tsession.query(Cassette_Connector)\
            .filter(and_(Cassette_Connector.creator == table_info['Connectors']['5'][0],
                         Cassette_Connector.creator_entry_number == table_info['Connectors']['5'][1])
                    )

        for left_query in left_connector:
            left_connector_overhang = left_query.overhang
        for right_query in right_connector:
            right_connector_overhang = right_query.overhang

        new_cassette_plasmid_entry = {'creator' : user_input.creator,
                                      'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                      'left_connector_part' : '1',
                                      'left_overhang' : left_connector_overhang,
                                      'right_connector_part' : '5',
                                      'right_overhang' : right_connector_overhang
                                      }

        Cassette_Plasmid.add(self.tsession, new_cassette_plasmid_entry, silent=False)

        for constituent_part in table_info['Part list ID']:
            new_cassette_assembly_entry = {'Cassette_creator' : user_input.creator,
                                           'Cassette_creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                           'Part_number' : constituent_part[0],
                                           'Part_creator' : constituent_part[1],
                                           'Part_creator_entry_number' : constituent_part[2]
                                           }
            Cassette_Assembly.add(self.tsession, new_cassette_assembly_entry, silent=False)

        self.add_features(table_info, current_plasmid_entry)
        self.upload_file(user_input, table_info, current_plasmid_entry)

        self.tsession.commit()

    def add_features(self, table_info, current_plasmid_entry):
        features_query = self.tsession.query(Feature, Feature_Type).filter(
            Feature.Feature_type == Feature_Type.Feature_type)
        possible_features = [
            (record.Feature_name, record.Feature_type, record.Feature_sequence, color.color, record.description) for
            record, color in features_query]
        for site in possible_features:
            target = site[2].upper()
            if table_info['Complete Assembly'].find(target) != -1:
                input_dict = {'creator' : current_plasmid_entry.creator,
                              'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                              'feature_name' : site[0]
                              }
                Plasmid_Feature.add(self.tsession, input_dict)

    def upload_file(self, user_input, table_info, current_plasmid_entry):
        buffy = self.generate_ape_file(user_input.database_ID, table_info['Complete Assembly'], table_info['Complete Description'])

        new_plasmid_file_entry = {'creator' : current_plasmid_entry.creator,
                                  'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                  'file_name' : '%s.gb' % user_input.UID,
                                  'file_type' : 'GenBank file',
                                  'Description' : current_plasmid_entry.description,
                                  'File' : buffy
                                  }

        Plasmid_File.add(self.tsession, new_plasmid_file_entry, silent=False)

    def generate_ape_file(self, plasmid_name, compete_assembly, complete_description, mutations = None, mutant_feature_tuple = None):
        features_list = []
        features_query = self.tsession.query(Feature, Feature_Type).filter(Feature.Feature_type == Feature_Type.Feature_type)
        possible_features = [(record.Feature_name, record.Feature_type, record.Feature_sequence, color.color, record.description) for record, color in features_query]

        # Special cases for mutations
        if mutant_feature_tuple != None:
            possible_features.append(mutant_feature_tuple)

        for site in possible_features:
            target = site[2].upper()
            #Find forward sequences from features database
            if compete_assembly.find(target) != -1:
                # Treat type II restriction enzymes as special cases
                if site[1] == 'Restrxn Type II':
                    # Forward Restriction Site
                    features_list.append(
                        SeqFeature(
                            CompoundLocation(
                                [
                                    FeatureLocation(compete_assembly.find(target),compete_assembly.find(target) + len(target)),
                                    FeatureLocation(compete_assembly.find(target) + 7, compete_assembly.find(target) + 11)
                                ]
                            ),
                            type = site[1],
                            strand = 1,
                            qualifiers = {"label": site[0],
                                          "ApEinfo_label": site[0],
                                          "ApEinfo_fwdcolor": site[3],
                                          "ApEinfo_revcolor":site[3],
                                          "ApEinfo_graphicformat":"arrow_data {{0 1 2 0 0 -1} {} 0}"
                                          }
                        )
                    )
                    # Reverse Restriction Site
                    features_list.append(
                            SeqFeature(
                                CompoundLocation(
                                    [
                                        FeatureLocation(compete_assembly.find(self.reverse_complement(target)),compete_assembly.find(self.reverse_complement(target)) + len(self.reverse_complement(target))),
                                        FeatureLocation(compete_assembly.find(self.reverse_complement(target)) - 5, compete_assembly.find(self.reverse_complement(target)) - 1)
                                    ]
                                ),
                                type = site[1],
                                strand = -1,
                                qualifiers = {"label": site[0],
                                              "ApEinfo_label": site[0],
                                              "ApEinfo_fwdcolor": site[3],
                                              "ApEinfo_revcolor":site[3],
                                              "ApEinfo_graphicformat":"arrow_data {{0 1 2 0 0 -1} {} 0}"
                                              }
                            )
                        )
                else:
                    features_list.append(
                        SeqFeature(
                            FeatureLocation(compete_assembly.find(target), compete_assembly.find(target) + len(target)),
                            type=site[1],
                            strand=1,
                            qualifiers={"label": site[0],
                                        "ApEinfo_label": site[0],
                                        "ApEinfo_fwdcolor": site[3],
                                        "ApEinfo_revcolor": site[3],
                                        "ApEinfo_graphicformat": "arrow_data {{0 1 2 0 0 -1} {} 0}"
                                        }
                        )
                    )

            if compete_assembly.find(self.reverse_complement(target)) != -1:
                if site[1] == 'Restrxn Type II':
                    pass
                else:
                    features_list.append(
                        SeqFeature(
                            FeatureLocation(compete_assembly.find(self.reverse_complement(target)),
                                            compete_assembly.find(self.reverse_complement(target)) + len(self.reverse_complement(target))),
                            type=site[1],
                            strand=-1,
                            qualifiers={"label": site[1],
                                        "ApEinfo_label": site[0],
                                        "ApEinfo_fwdcolor": site[3],
                                        "ApEinfo_revcolor": site[3],
                                        "ApEinfo_graphicformat": "arrow_data {{0 1 2 0 0 -1} {} 0}"
                                        }
                        )
                    )


        if mutations != None:
            for mutation in mutations:
                features_list.append(mutation)

        sequence = SeqRecord( Seq(compete_assembly,
                                  IUPAC.unambiguous_dna),
                              id=plasmid_name,
                              name=plasmid_name,
                              description=complete_description,
                              features = features_list
                              )

        buffy = StringIO.StringIO()
        SeqIO.write(sequence, buffy, 'gb')

        return buffy.getvalue()

    def mutant_check(self, temp_sequence):
        sequence_pass = True
        if 'GGTCTC' in temp_sequence \
                or 'GAGACC' in temp_sequence \
                or 'CGTCTC' in temp_sequence \
                or 'GAGACG' in temp_sequence:
            sequence_pass = False
        return sequence_pass

    def generate_mutant_sequence(self, tsession, user_mutations):
        # E. coli codon table: http://www.sci.sdsu.edu/~smaloy/MicrobialGenetics/topics/in-vitro-genetics/codon-usage.html

        codon_table = {'A' : ['GCG', 'GCC', 'GCA', 'GCT'],
                       'C' : ['TGC', 'TGT'],
                       'D' : ['GAT', 'GAC'],
                       'E' : ['GAA', 'GAG'],
                       'F' : ['TTT', 'TTC'],
                       'G' : ['GGC', 'GGT', 'GGG', 'GGA'],
                       'H' : ['CAT', 'CAC'],
                       'I' : ['ATT', 'ATC', 'ATA'],
                       'K' : ['AAA', 'AAG'],
                       'L' : ['CTG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA'],
                       'M' : ['AUG'],
                       'N' : ['AAC', 'AAT'],
                       'P' : ['CCG', 'CCA', 'CCT', 'CCC'],
                       'Q' : ['CAG', 'CAA'],
                       'R' : ['CGT', 'CGC', 'CGG', 'CGA', 'AGA', 'AGG'],
                       'S' : ['AGC', 'TCT', 'TCC', 'TCG', 'AGT', 'TCA'],
                       'T' : ['ACC', 'ACA', 'ACG', 'ACU'],
                       'V' : ['GTG', 'GTT' ,'GTC', 'GTA'],
                       'W' : ['TGG'],
                       'Y' : ['TAT', 'TAC']
                       }

        feature_query = tsession.query(Plasmid_Feature, Feature).filter(Plasmid_Feature.ID == user_mutations.Plasmid_Feature_ID).filter(Plasmid_Feature.feature_name == Feature.Feature_name )

        for plasmid_feature, feature in feature_query:
            Actual_WT_feature_sequence = feature.Feature_sequence
            WT_feature_sequence = feature.Feature_sequence
            # print plasmid_feature.ID
            # print plasmid_feature.creator
            # print plasmid_feature.creator_entry_number
            # print feature.Feature_name
            # print feature.Feature_type
            # print feature.Feature_sequence


        for point_mutant in user_mutations.mutation_list:

            assert WT_feature_sequence[3 * point_mutant[1] - 3 : -(len(WT_feature_sequence) - 3 * point_mutant[1])] in codon_table[point_mutant[0]], \
                'User input target codon and sequence codon identities do not match!'

            mutant_codon_choices = codon_table[point_mutant[2]]

            for codon in mutant_codon_choices:
                temp_sequence = WT_feature_sequence[ : 3 * point_mutant[1] - 3 ] + codon + WT_feature_sequence[ - ( len(WT_feature_sequence) - ( 3 * point_mutant[1] ) ) : ]
                if self.mutant_check(temp_sequence):
                    break
                else:
                    continue
            assert self.mutant_check(temp_sequence), "No suitable codons were found for %s%s%s! You're going to have to do this one by hand :P" %(point_mutant[0], point_mutant[1], point_mutant[2])

            WT_feature_sequence = temp_sequence.upper()

        sequence_query = tsession.query(Plasmid_Feature, Feature, Feature_Type, Plasmid)\
            .filter(Plasmid_Feature.ID == user_mutations.Plasmid_Feature_ID)\
            .filter(Plasmid.creator == Plasmid_Feature.creator)\
            .filter(Plasmid.creator_entry_number == Plasmid_Feature.creator_entry_number)\
            .filter(Plasmid_Feature.feature_name == Feature.Feature_name)\
            .filter(Feature.Feature_type == Feature_Type.Feature_type)

        for plasmid_feature, feature, feature_type, plasmid in sequence_query:
            database_ID = 'p%s%s' %(plasmid.creator, ('0000' + str(plasmid.creator_entry_number))[-4:])
            WT_plasmid_sequence = plasmid.sequence
            WT_feature_name = plasmid_feature.feature_name
            WT_feature_color = feature_type.color
            WT_feature_type = feature_type.Feature_type
            WT_feature_description = feature.description

            # Tuple for annotating mutant version of feature
            mutant_feature_tuple = ('Mutant %s' % WT_feature_name,
                                    WT_feature_type,
                                    WT_feature_sequence,
                                    WT_feature_color,
                                    WT_feature_description
                                    )

        print WT_feature_sequence

        # Generate plasmid sequence with mutant CDS
        CDS_index_start = WT_plasmid_sequence.find(Actual_WT_feature_sequence)
        Mutant_plasmid_sequence = WT_plasmid_sequence[ : CDS_index_start ] + WT_feature_sequence + WT_plasmid_sequence[ -( len(WT_plasmid_sequence) - (CDS_index_start + len(WT_feature_sequence)) ) : ]

        print WT_plasmid_sequence[ : CDS_index_start ]
        print '\n'
        print WT_feature_sequence
        print '\n'
        print WT_plasmid_sequence[ -( len(WT_plasmid_sequence) - (CDS_index_start + len(WT_feature_sequence)) ) : ]


        # Annotate mutant codons
        Mutant_codon_list = []

        for mutation in user_mutations.mutation_list:

            Mutant_codon_list.append(
                SeqFeature(
                    FeatureLocation(CDS_index_start + mutation[1] * 3 - 3, CDS_index_start + mutation[1] * 3),
                    type='Mutation',
                    strand=1,
                    qualifiers={"label": '%s%s%s' %(mutation[0], mutation[1], mutation[2]),
                                "ApEinfo_label": '%s%s%s' %(mutation[0], mutation[1], mutation[2]),
                                "ApEinfo_fwdcolor": '#506380',
                                "ApEinfo_revcolor": '#506380',
                                "ApEinfo_graphicformat": "arrow_data {{0 1 2 0 0 -1} {} 0}"
                                }
                )
            )

        mutant_genbank_file = self.generate_ape_file((database_ID + '_Mutant'),
                                                     Mutant_plasmid_sequence,
                                                     WT_feature_description,
                                                     mutations=Mutant_codon_list,
                                                     mutant_feature_tuple=mutant_feature_tuple)

        with open('Mutant_GenBank_Test.gb', 'w+b') as file:
            file.write(mutant_genbank_file)




