import sys

sys.path.insert(0, '..')

from db.interface import DatabaseInterface
from db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature


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
    * Upload assembled plasmids to the database

    ###########################
    # Planned Funtionality:   #
    ###########################

    * Assemble and submit part plasmids to the database

    '''
    def __init__(self):
        # Create up the database session
        self.dbi = DatabaseInterface()
        self.tsession = self.dbi.get_session()

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

        :param part_plasmid_list: A list of part plasmids that the user wants to assemble into a cassette plasmid
        :param assembly_type: Cassette or Multicassette
        :return: complete plasmid assembly as string
        '''
        # print plasmid_list

        def fetch_cassette_parts():
            # Get sequences and part numbers for listed part plasmids
            part_plasmids_query = self.tsession.query(Plasmid, Part_Plasmid, Part_Plasmid_Part, Part_Type)\
                .filter(Part_Plasmid.creator == Plasmid.creator)\
                .filter(Part_Plasmid.creator_entry_number == Plasmid.creator_entry_number)\
                .filter(Part_Plasmid.creator == Part_Plasmid_Part.creator)\
                .filter(Part_Plasmid.creator_entry_number == Part_Plasmid_Part.creator_entry_number)\
                .filter(Plasmid.plasmid_name.in_(user_input.input_plasmids)) \
                .filter(Part_Plasmid_Part.part_number == Part_Type.part_number)

            table_info = {}
            table_info['Part list'] = []
            table_info['Part list ID'] = []
            table_info['Description'] = []
            table_info['Connectors'] = {}
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

        # Start with Part one (cassette assembly) or backbone (multicassette assembly) so that there is consistency in how the plasmids are assembled
        if user_input.assembly_type == 'Cassette':
            table_info = fetch_cassette_parts()
            for part in table_info['Part list']:
                if part[:4].upper() == 'CCCT':
                    intermediate = part.upper()
                    table_info['Part list'].remove(part)
                    break

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

        table_info['Complete Assembly'] = intermediate[:-4].upper()
        table_info['Complete Description'] = '\n'.join(table_info['Description'])

        return table_info

    def plasmid_checks(self, table_info, user_input):
        '''
        Verifies that a plasmid conforms to modular cloning sequence standards before entry into database

        :return:
        '''
        if user_input.assembly_type == 'Part':
            pass

        if user_input.assembly_type == 'Cassette':
            assert table_info['Complete Assembly'].count('CGTCTC') == 1, \
                'There is more than one forward %s site in %s! This plasmid contains %s forward %s sites.' % ('BsmBI', user_input.UID, table_info['Complete Assembly'].count('CGTCTC'), 'BsmBI')
            assert table_info['Complete Assembly'].count('GAGACG') == 1, \
                'There is more than one reverse %s site in %s! This plasmid contains %s reverse %s sites.' % ('BsmBI', user_input.UID, table_info['Complete Assembly'].count('GAGACG'), 'BsmBI')

        if user_input.assembly_type == 'Multicassette':
            assert table_info['Complete Assembly'].count('GGTCTC') == 1, \
                'There is more than one forward %s site in %s! This plasmid contains %s forward %s sites.' % (
                'BsaI', user_input.UID, table_info['Complete Assembly'].count('GGTCTC'), 'BsaI')
            assert table_info['Complete Assembly'].count('GAGACC') == 1, \
                'There is more than one reverse %s site in %s! This plasmid contains %s reverse %s sites.' % (
                'BsaI', user_input.UID, table_info['Complete Assembly'].count('GAGACC'), 'BsaI')

    def add_cassette_plasmid_to_db(self, user_input, table_info):
        new_plasmid_entry = Plasmid(creator = user_input.creator,
                                    plasmid_name = user_input.UID,
                                    plasmid_type = user_input.assembly_type,
                                    location = user_input.location,
                                    description = table_info['Complete Description'],
                                    sequence = table_info['Complete Assembly'],
                                    status = 'designed'
                                    )

        self.tsession.add(new_plasmid_entry)
        self.tsession.flush()

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

        new_cassette_plasmid_entry = Cassette_Plasmid(creator = user_input.creator,
                                                      creator_entry_number = new_plasmid_entry.creator_entry_number,
                                                      left_connector_part = '1',
                                                      left_overhang = left_connector_overhang,
                                                      right_connector_part = '5',
                                                      right_overhang = right_connector_overhang
                                                      )

        self.tsession.add(new_cassette_plasmid_entry)
        self.tsession.flush()

        for constituent_part in table_info['Part list ID']:
            new_cassette_assembly_entry = Cassette_Assembly(Cassette_creator = user_input.creator,
                                                            Cassette_creator_entry_number = new_plasmid_entry.creator_entry_number,
                                                            Part_number = constituent_part[0],
                                                            Part_creator = constituent_part[1],
                                                            Part_creator_entry_number = constituent_part[2]
                                                            )
            self.tsession.add(new_cassette_assembly_entry)

        self.tsession.flush()
        # self.tsession.commit()

    def generate_ape_file(self, user_input, table_info):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        from Bio.Alphabet import IUPAC

        features_list = []
        features_query = self.tsession.query(Feature, Feature_Type).filter(Feature.Feature_type == Feature_Type.Feature_type)
        possible_features = [(record.Feature_name, record.Feature_type, record.Feature_sequence, color.color, record.description) for record, color in features_query]

        for site in possible_features:
            target = site[2].upper()
            #Find forward sequences from features database
            if table_info['Complete Assembly'].find(target) != -1:
                # Treat type II restriction enzymes as special cases
                if site[1] == 'Restrxn Type II':
                    # Forward Restriction Site
                    features_list.append(
                        SeqFeature(
                            CompoundLocation(
                                [
                                    FeatureLocation(table_info['Complete Assembly'].find(target),table_info['Complete Assembly'].find(target) + len(target)),
                                    FeatureLocation(table_info['Complete Assembly'].find(target) + 7, table_info['Complete Assembly'].find(target) + 11)
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
                                        FeatureLocation(table_info['Complete Assembly'].find(self.reverse_complement(target)),table_info['Complete Assembly'].find(self.reverse_complement(target)) + len(self.reverse_complement(target))),
                                        FeatureLocation(table_info['Complete Assembly'].find(self.reverse_complement(target)) - 5, table_info['Complete Assembly'].find(self.reverse_complement(target)) - 1)
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
                            FeatureLocation(table_info['Complete Assembly'].find(target), table_info['Complete Assembly'].find(target) + len(target)),
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

            if table_info['Complete Assembly'].find(self.reverse_complement(target)) != -1:
                if site[1] == 'Restrxn Type II':
                    pass
                else:
                    features_list.append(
                        SeqFeature(
                            FeatureLocation(table_info['Complete Assembly'].find(self.reverse_complement(target)),
                                            table_info['Complete Assembly'].find(self.reverse_complement(target)) + len(self.reverse_complement(target))),
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


        # features_list.append(SeqFeature(FeatureLocation(100,200), type = 'asdf', strand = 1, qualifiers={'label':'qwer'}))

        sequence = SeqRecord( Seq(table_info['Complete Assembly'],
                                  IUPAC.unambiguous_dna),
                              id=user_input.UID,
                              name=user_input.UID,
                              description=table_info['Complete Description'],
                              features = features_list
                              )

        with open('%s-TEST.gb' % user_input.UID, 'w') as final_ape_output:
            SeqIO.write(sequence, final_ape_output, 'gb')





