import sys
import StringIO
import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC
from sqlalchemy import and_, or_

sys.path.insert(0, '..')

try:
    from db.interface import DatabaseInterface
    from db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature, Plasmid_Feature, Plasmid_File, CDS_Mutant, CDS_Mutant_Constituent, Multicassette_Assembly, Multicassette_Plasmid, Publication_Plasmid, Other_Plasmid, Plasmid_Feature_Design, Design_Plasmid, Plasmid_Design_Ancestors, Plasmid_Design_Features
    from db.procedures import call_procedure
except:
    # nasty hack since we are not packaging things up properly yet for external use (e.g. the website)
    from kprimers.db.interface import DatabaseInterface
    from kprimers.db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature, Plasmid_Feature, Plasmid_File, CDS_Mutant, CDS_Mutant_Constituent, Multicassette_Assembly, Multicassette_Plasmid, Publication_Plasmid, Other_Plasmid, Plasmid_Feature_Design, Design_Plasmid, Plasmid_Design_Ancestors, Plasmid_Design_Features
    from kprimers.db.procedures import call_procedure
class Plasmid_Exception(Exception): pass

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
    * Storing information for mutations and libraries
    * Generate mutant sequences for plasmid features and associated .ape file
    * Multicassette assembly and submission

    ###########################
    # Planned Funtionality:   #
    ###########################

    * Annotating working and archive locations for all plasmids

    '''
    def __init__(self, tsession = None, username = None):
        # Create up the database session
        self.dbi = None
        try:
            self.dbi = DatabaseInterface()
        except:
            if not tsession:
                raise
        self.tsession = tsession or self.dbi.get_session()

        # Get User ID
        self.username = username
        if not self.username:
            self.username = u'' + os.getlogin() # The webserver cannot call this function so it must instead pass username in as an argument
        self.user_ID = self.tsession.query(Users).filter(Users.lab_username == self.username).one().ID


    def __del__(self):
        if self.tsession:
            self.tsession.close()


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

        # Create the database session
        dbi = self.dbi
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

    def golden_gate_assembly(self, input_sequences, assembly_type, part_type=None):
        '''
        This function will take user input of a list of plasmids and attempt to perform a golden gate assembly.

        :return:
        '''

        table_info = {}
        table_info['Part list'] = []
        table_info['Part list ID'] = []
        table_info['Description'] = []
        table_info['Connectors'] = {}

        if assembly_type.lower() == 'part':
            # Checks that part 3s have len % 3 = 0 and do not contain a starting MET or stop codon anywhere
            if '3' in part_type:
                if len(input_sequences[0]) % 3 != 0:
                    raise Plasmid_Exception('Your input coding sequence contains incomplete codons!!')

            if any(input_type == part_type for input_type in ['3', '3a']):
                # Checks that input sequence does not start with MET
                # Update: removed constraint on MET, added warning instead
                if input_sequences[0][:3].upper() == 'ATG':
                    print ('WARNING: The start codon is already in the part plasmid sequence, you do not need to include a MET start codon.')

                # Checks that stop codons are not in the input coding sequence
                if any(codon in [input_sequences[0][i:i + 3].upper() for i in range(0, len(input_sequences[0]), 3)] for codon in ['TAA', 'TAG', 'TGA']):
                    raise Plasmid_Exception('There are stop codons in your coding sequence! Please remove them!')

            part_entry_vector = self.tsession.query(Plasmid).filter(Plasmid.creator == u'JL').filter(Plasmid.creator_entry_number == 2).one()
            sequence_upper = self.add_part_arms(input_sequences[0], part_type)
            table_info['Part list'].append(sequence_upper[ sequence_upper.find('CGTCTC') + 7 : sequence_upper.find('GAGACG') - 1])
            intermediate = part_entry_vector.sequence[ : (part_entry_vector.sequence.upper().find('GAGACG') - 1) ]

        # Start with Part one (cassette assembly) so that there is consistency in how the plasmids are assembled
        if assembly_type.lower() == 'cassette':
            table_info = self.fetch_cassette_parts(input_sequences, table_info)
            for part in table_info['Part list']:
                if part[:4].upper() == 'CCCT':
                    intermediate = part.upper()
                    table_info['Part list'].remove(part)
                    break

        # Start with Part one (cassette assembly) so that there is consistency in how the plasmids are assembled
        if assembly_type.lower() == 'multicassette':
            table_info = self.fetch_mulitcassette_cassettes(input_sequences, table_info)

            for part in table_info['Part list']:
                if part[:4].upper() == 'CTGA':
                    intermediate = part.upper()
                    table_info['Part list'].remove(part)
                    break

        assembly_step = 0
        assembly_steps = len(table_info['Part list'])

        while assembly_step < assembly_steps:
            for part in table_info['Part list']:
                # If part 5' and intermediate 3' match...
                if part.upper()[:4] == intermediate.upper()[-4:]:
                    intermediate = intermediate.upper()[:-4] + part.upper()
            assembly_step += 1

        if intermediate[-4:] != intermediate[:4]:
            raise Plasmid_Exception('Incomplete assembly! This assembly does not produce a circular plasmid! :(')

        table_info['Complete Assembly'] = intermediate[4:].upper()
        table_info['Complete Description'] = ' | '.join(table_info['Description'])

        return table_info


    def add_part_arms(self, sequence, part_type):
        # Some of the overhangs aren't added yet since I don't have access to APE at the moment...
        arm_dict = {'1': ['', ''],
                    '2': ['gcatCGTCTCaAGCAGGTCTCAAACG', 'TATGtGAGACCtGAGACGgcat'],
                    # '2a' : ['', ''],
                    #  '2b' : ['', ''],
                    '3': ['gcatCGTCTCaAGCAGGTCTCaTATG', 'taaATCCtGAGACCtGAGACGgcat'],
                    '3a': ['gcatCGTCTCaAGCAGGTCTCaTATG', 'GGTAGCGGCAGCGGCAGTTCTtGAGACCtGAGACGgcat'],
                    '3b': ['gcatCGTCTCaAGCAGGTCTCATTCT', 'taaATCCtGAGACCtGAGACGgcat'],
                    '4': ['gcatCGTCTCaAGCAGGTCTCaATCC', 'GCTGtGAGACCtGAGACGgcat'],
                    # '4a' : ['gcatCGTCTCaAGCAGGTCTCaATCC', ''],
                    #  '4b' : ['', 'GCTGtGAGACCtGAGACGgcat'],
                    '5': ['', ''],
                    '6': ['gcatCGTCTCaAGCAGGTCTCA', 'TGAGACCtGAGACGgcat'],
                    '7': ['', '']
                    }

        try:
            part_sequence = arm_dict[part_type][0] + sequence + arm_dict[part_type][1]
            if part_sequence.upper().find('CGTCTC') == -1 or part_sequence.upper().find('GAGACG') == -1:
                raise Plasmid_Exception('BsmBI sites not found! Something fishy is going on...')

        except:
            raise Plasmid_Exception('%s part types are not supported yet!' % part_type)

        return part_sequence.upper()


    def fetch_cassette_parts(self, input_sequences, table_info=None):
        # Get sequences and part numbers for listed part plasmids

        part_IDs = [and_(Plasmid.creator == part_creator, Plasmid.creator_entry_number == part_creator_entry_number) for (part_creator, part_creator_entry_number) in input_sequences]

        part_plasmids_query = self.tsession.query(Plasmid, Part_Plasmid, Part_Plasmid_Part, Part_Type) \
            .filter(and_(Part_Plasmid.creator == Plasmid.creator,
                         Part_Plasmid.creator_entry_number == Plasmid.creator_entry_number,
                         Part_Plasmid.creator == Part_Plasmid_Part.creator,
                         Part_Plasmid.creator_entry_number == Part_Plasmid_Part.creator_entry_number,
                         Part_Plasmid_Part.part_number == Part_Type.part_number,
                         or_(*part_IDs)
                         )
                    )

        restriction_enzyme = 'BsaI'
        site_F = 'GGTCTC'
        site_R = 'GAGACC'

        already_fetched_cassettes = []
        part_sequences = []

        for plasmid, part_plasmid, part_plasmid_part, part_type in part_plasmids_query:

            print "%s\t%s\t%s\t%s" % (
            plasmid.plasmid_name, plasmid.creator, plasmid.creator_entry_number, part_plasmid_part.part_number)

            sequence_upper = plasmid.sequence.upper()

            if sequence_upper.count(site_F) != 1:
                raise Plasmid_Exception('There is more than one forward %s site in %s! %s' % (restriction_enzyme, plasmid.plasmid_name or plasmid.get_id(), plasmid.sequence.count(site_F)))
            if sequence_upper.count(site_R) != 1:
                raise Plasmid_Exception('There is more than one reverse %s site in %s!' % (restriction_enzyme, plasmid.plasmid_name or plasmid.get_id()))

            left_overhang = part_type.overhang_5
            right_overhang = part_type.overhang_3

            # check if a part plasmid has multiple part plasmid parts - will select the correct 5' and 3' overhangs
            if part_plasmids_query.filter( Plasmid.creator == plasmid.creator).filter(Plasmid.creator_entry_number == plasmid.creator_entry_number).count() > 1:
                left_overhangs = set()
                right_overhangs = set()

                for p, pp, ppp, pt in part_plasmids_query.filter(and_(Plasmid.creator == plasmid.creator, Plasmid.creator_entry_number == plasmid.creator_entry_number)):
                    left_overhangs.add(pt.overhang_5)
                    right_overhangs.add(pt.overhang_3)

                leftoverhangs = left_overhangs.symmetric_difference(right_overhangs) # Leftover overhangs? Get it?

                for p, pp, ppp, pt in part_plasmids_query.filter(and_(Plasmid.creator == plasmid.creator, Plasmid.creator_entry_number == plasmid.creator_entry_number)):
                    if pt.overhang_5 in leftoverhangs:
                        left_overhang = pt.overhang_5
                    if pt.overhang_3 in leftoverhangs:
                        right_overhang = pt.overhang_3

            site_F_position = sequence_upper.find(left_overhang, sequence_upper.find(site_F))
            site_R_position = sequence_upper.find(right_overhang, sequence_upper.find(site_R) - 10) + 4

            # Circular permutation so that forward cut site is always upstream of reverse cut site in linear sequence
            if site_F_position > site_R_position:
                sequence_upper = sequence_upper[site_R_position + 8:] + sequence_upper[:site_R_position + 8]

            if table_info:
                if (plasmid.creator, plasmid.creator_entry_number) not in already_fetched_cassettes:
                    table_info['Part list'].append(sequence_upper[
                                                   sequence_upper.find(left_overhang, sequence_upper.find(site_F)):
                                                   sequence_upper.find(right_overhang,
                                                                       sequence_upper.find(site_R) - 10) + 4
                                                   ])
                    table_info['Description'].append(plasmid.description)
                    table_info['Part list ID'].append(
                        (part_plasmid_part.part_number, plasmid.creator, plasmid.creator_entry_number))
                    if part_plasmid_part.part_number == '1' or part_plasmid_part.part_number == '5':
                        table_info['Connectors'][part_plasmid_part.part_number] = [plasmid.creator,
                                                                                   plasmid.creator_entry_number]
                    already_fetched_cassettes.append((plasmid.creator, plasmid.creator_entry_number))

            new_part_sequence = sequence_upper[sequence_upper.find(left_overhang, sequence_upper.find(site_F)):sequence_upper.find(right_overhang,sequence_upper.find(site_R) - 10) + 4]
            if new_part_sequence not in part_sequences:
                part_sequences.append(new_part_sequence)

        if table_info:
            return table_info
        else:
            return part_sequences


    def fetch_mulitcassette_cassettes(self, input_sequences, table_info=None):
        # Get sequences for cassettes in input_sequences
        restriction_enzyme = 'BsmBI'
        site_F = 'CGTCTC'
        site_R = 'GAGACG'

        cassettee_IDs = [and_(Plasmid.creator == cassette_creator, Plasmid.creator_entry_number == cassette_creator_entry_number) for
                    (cassette_creator, cassette_creator_entry_number) in input_sequences]

        cassette_plasmids_query = self.tsession.query(Plasmid, Cassette_Plasmid) \
            .filter(and_(Cassette_Plasmid.creator == Plasmid.creator,
                         Cassette_Plasmid.creator_entry_number == Plasmid.creator_entry_number,
                         or_(*cassettee_IDs)
                         )
                    )

        part_sequences = []

        for plasmid, cassette_plasmid in cassette_plasmids_query:

            print "%s\t%s\t%s\t%s" % (
                cassette_plasmid.creator, cassette_plasmid.creator_entry_number, cassette_plasmid.left_overhang, cassette_plasmid.right_overhang)

            sequence_upper = plasmid.sequence.upper()

            if sequence_upper.count(site_F) != 1:
                raise Plasmid_Exception('There is more than one forward %s site in %s! %s' % (restriction_enzyme, plasmid.plasmid_name or plasmid.get_id(), plasmid.sequence.count(site_F)))
            if sequence_upper.count(site_R) != 1:
                raise Plasmid_Exception('There is more than one reverse %s site in %s!' % (restriction_enzyme, plasmid.plasmid_name or plasmid.get_id()))

            left_overhang = cassette_plasmid.left_overhang
            right_overhang = cassette_plasmid.right_overhang

            site_F_position = sequence_upper.find(left_overhang, sequence_upper.find(site_F))
            site_R_position = sequence_upper.find(right_overhang, sequence_upper.find(site_R) - 10) + 4

            # Circular permutation so that forward cut site is always upstream of reverse cut site in linear sequence
            if site_F_position > site_R_position:
                sequence_upper = sequence_upper[site_R_position + 8:] + sequence_upper[:site_R_position + 8]

            if table_info:
                table_info['Part list'].append(sequence_upper[
                                               sequence_upper.find(left_overhang, sequence_upper.find(site_F)):
                                               sequence_upper.find(right_overhang,
                                                                   sequence_upper.find(site_R) - 10) + 4
                                               ])
                table_info['Description'].append(plasmid.description)
                table_info['Part list ID'].append((plasmid.creator, plasmid.creator_entry_number))
            part_sequences.append(sequence_upper[sequence_upper.find(left_overhang, sequence_upper.find(site_F)):sequence_upper.find(right_overhang, sequence_upper.find(site_R) - 10) + 4])

        if table_info:
            return table_info
        else:
            return part_sequences


    def plasmid_checks(self, input_dict, plasmid_or_assembly_type, part_type = None):
        '''
        Verifies that a plasmid conforms to modular cloning sequence standards before entry into database

        # Part Plasmids #
          * Part 1 and Part 5 should contain one reverse BsmBI site
          * If not a Part 1 or Part 5, sequence should not contain BsmBI sites
          *
        :return:
        '''

        #todo: write checks for other plasmids? If that exists?

        if plasmid_or_assembly_type.lower() == 'part':
            # If part plasmid contains multiple part_plasmid_parts
            # if part_type == '1' or part_type == '5':
            if '1' in re.split(",| ", part_type) or '5' in re.split(",| ", part_type):
                if input_dict['sequence'].count('GAGACG') != 1:
                    raise Plasmid_Exception('There should be only one reverse %s site in a Connector Part! Your part %s submission contains %s reverse %s sites.' % (
                    'BsmBI', part_type, input_dict['sequence'].count('GAGACG'), 'BsmBI'))
            else:
                if input_dict['sequence'].count('GAGACG') != 0:
                    raise Plasmid_Exception('There should not be any reverse %s sites in this part type! Your submission %s contains %s reverse %s site(s).' % (
                        'BsmBI', input_dict['plasmid_name'], input_dict['sequence'].count('GAGACG'), 'BsmBI'))
            if input_dict['sequence'].count('CGTCTC') != 0:
                raise Plasmid_Exception('There should not be any forward %s sites in this part type! This plasmid contains %s forward %s sites.' % (
                    'BsmBI', input_dict['sequence'].count('CGTCTC'), 'BsmBI'))
            if input_dict['sequence'].count('GGTCTC') != 1:
                raise Plasmid_Exception('There is more than one forward %s site in %s! Your submission contains %s forward %s sites.' % (
                    'BsaI', input_dict['plasmid_name'], input_dict['sequence'].count('GGTCTC'), 'BsaI'))
            if input_dict['sequence'].count('GAGACC') != 1:
                raise Plasmid_Exception('There is more than one reverse %s site in %s! Your submission contains %s reverse %s sites.' % (
                    'BsaI', input_dict['plasmid_name'], input_dict['sequence'].count('GAGACC'), 'BsaI'))

        elif plasmid_or_assembly_type.lower() == 'cassette':
            if input_dict['sequence'].count('GGTCTC') != 0:
                raise Plasmid_Exception('There are BsaI sites in your completed cassette assembly! This thing is going to get wrecked during the final digestion step of the Golden Gate Assembly!')
            if input_dict['sequence'].count('GAGACC') != 0:
                raise Plasmid_Exception('There are BsaI sites in your completed cassette assembly! This thing is going to get wrecked during the final digestion step of the Golden Gate Assembly!')

            #todo: Uncomment this entire section and make the correct checks... comments below were made assuming the ConLS connector was correct... it was not.
            # I need to implement more complicated logic when checking for cassette connector sites in cassette assemblies
            # It looks like LS contains a forward BsmBI site and all other connectors have a reverse BsmBI site...
            # Not going to deal with that right now...

            # if input_dict['sequence'].count('CGTCTC') != 0:
            #     raise Plasmid_Exception('Your submission contains %s forward %s sites. There should be none!' % ( input_dict['sequence'].count('CGTCTC'), 'BsmBI'))
            # if input_dict['sequence'].count('GAGACG') != 2:
            #     raise Plasmid_Exception('There should only be two reverse %s sites in a cassette assembly. Your submission contains %s reverse %s sites.' % ('BsmBI', input_dict['sequence'].count('GAGACG'), 'BsmBI'))

        elif plasmid_or_assembly_type.lower() == 'multicassette':
            # Ehhh... I'll leave this for now. There really aren't any restrictions for a multicasssette plasmid that I can think of for now
            if input_dict['sequence'].count('CGTCTC') != 0:
                raise Plasmid_Exception('There is a forward BsmBI site in your assembled plasmid. Get rid of it!!!')
            if input_dict['sequence'].count('GAGACG') != 0:
                raise Plasmid_Exception('There is a reverse BsmBI site in your assembled plasmid. Get rid of it!!!')

        else:
            # todo: 'other' plasmid type designation will need to be updated, changed to 'standard' for the time being
            # This forces assembly_type to be passed correctly. Shane wrote a bug by not passing this argument in correctly.
            assert(plasmid_or_assembly_type.lower() == 'standard')


    def add_part_plasmid_to_db(self, input_dict, part_type, auto_commit = False):
        '''
        :param user_input:
        :param table_info:
        :param auto_commit:
        :return: The created Plasmid object.
        '''
        current_plasmid_entry = Plasmid.add(self.tsession, input_dict, silent=True)

        new_part_plasmid_entry = {'creator' : current_plasmid_entry.creator,
                                  'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                  'resistance' : u'CM'
                                  }
        Part_Plasmid.add(self.tsession, new_part_plasmid_entry, silent=True)

        for part in part_type.split(','):
            new_part_plasmid_part_entry = {'creator' : current_plasmid_entry.creator,
                                           'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                           'part_number' : part.strip()
                                           }
            Part_Plasmid_Part.add(self.tsession, new_part_plasmid_part_entry, silent = True)

        self.add_features(current_plasmid_entry)
        self.upload_file(current_plasmid_entry)

        if auto_commit:
            self.tsession.commit()

        return current_plasmid_entry


    def add_cassette_plasmid_to_db(self, input_dict, table_info, auto_commit = False):

        current_plasmid_entry = Plasmid.add(self.tsession, input_dict, silent=False)

        from sqlalchemy import and_
        left_connector = self.tsession.query(Cassette_Connector)\
            .filter(and_(Cassette_Connector.creator == table_info['Connectors']['1'][0],
                         Cassette_Connector.creator_entry_number == table_info['Connectors']['1'][1])
                    ).one()
        right_connector = self.tsession.query(Cassette_Connector)\
            .filter(and_(Cassette_Connector.creator == table_info['Connectors']['5'][0],
                         Cassette_Connector.creator_entry_number == table_info['Connectors']['5'][1])
                    ).one()

        left_connector_overhang = left_connector.overhang
        right_connector_overhang = right_connector.overhang

        new_cassette_plasmid_entry = {'creator' : current_plasmid_entry.creator,
                                      'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                      'left_connector_part' : '1',
                                      'left_overhang' : left_connector_overhang,
                                      'right_connector_part' : '5',
                                      'right_overhang' : right_connector_overhang
                                      }

        Cassette_Plasmid.add(self.tsession, new_cassette_plasmid_entry, silent=False)

        for constituent_part in table_info['Part list ID']:
            new_cassette_assembly_entry = {'Cassette_creator' : current_plasmid_entry.creator,
                                           'Cassette_creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                           'Part_number' : constituent_part[0],
                                           'Part_creator' : constituent_part[1],
                                           'Part_creator_entry_number' : constituent_part[2]
                                           }
            Cassette_Assembly.add(self.tsession, new_cassette_assembly_entry, silent=False)

        self.add_features(current_plasmid_entry)
        self.upload_file(current_plasmid_entry)

        if auto_commit:
            self.tsession.commit()

        return current_plasmid_entry


    def add_multicassette_plasmid_to_db(self, input_dict, table_info, auto_commit=False):

        current_plasmid_entry = Plasmid.add(self.tsession, input_dict, silent=False)
        mulitcassette_plasmid_added = False

        for constituent_part in table_info['Part list ID']:
            overhangs = self.tsession.query(Cassette_Plasmid).filter(and_(Cassette_Plasmid.creator == constituent_part[0], Cassette_Plasmid.creator_entry_number == constituent_part[1])).one()

            if overhangs.left_connector_part == 'BS' and overhangs.right_connector_part == 'BE':
                new_multicassette_plasmid_entry = {'creator': current_plasmid_entry.creator,
                                              'creator_entry_number': current_plasmid_entry.creator_entry_number,
                                              'Backbone_creator': overhangs.creator,
                                              'Backbone_creator_entry_number': overhangs.creator_entry_number,
                                              }

                Multicassette_Plasmid.add(self.tsession, new_multicassette_plasmid_entry, silent=False)
                mulitcassette_plasmid_added = True

            new_multicassette_assembly_entry = {'Multicassette_creator': current_plasmid_entry.creator,
                                                'Multicassette_creator_entry_number': current_plasmid_entry.creator_entry_number,
                                                'Left_connector' : overhangs.left_overhang,
                                                'Right_connector' : overhangs.right_overhang,
                                                'Cassette_creator': constituent_part[0],
                                                'Cassette_creator_entry_number': constituent_part[1]
                                                }

            Multicassette_Assembly.add(self.tsession, new_multicassette_assembly_entry, silent=False)

        if not mulitcassette_plasmid_added:
            raise

        self.add_features(current_plasmid_entry)
        self.upload_file(current_plasmid_entry)

        if auto_commit:
            self.tsession.commit()

        return current_plasmid_entry


    def add_other_plasmid_to_db(self, input_dict, input_files):
        new_plasmid = Plasmid.add(self.tsession, input_dict)

        # Add Plasmid File records
        plasmid_id_fields = dict(
            creator_entry_number=new_plasmid.creator_entry_number,
            creator=new_plasmid.creator,
        )
        for pfile in input_files:
            pfile.update(plasmid_id_fields)
            pf = Plasmid_File.add(self.tsession, pfile)

        self.add_features(new_plasmid)


    def add_features(self, current_plasmid_entry):
        # todo: add BsmBI and BsaI sites for MoClo plasmids only
        features_query = self.tsession.query(Feature, Feature_Type).filter(
            Feature.Feature_type == Feature_Type.Feature_type)
        possible_features = [
            (record.Feature_name, record.Feature_type, record.Feature_sequence, color.color, record.description) for
            record, color in features_query]
        for site in possible_features:
            target = site[2].upper()
            if current_plasmid_entry.sequence.upper().find(target) != -1 or current_plasmid_entry.sequence.upper().find(self.reverse_complement(target)) != -1:
                if site[1] == 'Restrxn Type II':
                    if current_plasmid_entry.plasmid_type in ['part', 'cassette', 'multicassette']:
                        input_dict = {'creator': current_plasmid_entry.creator,
                                      'creator_entry_number': current_plasmid_entry.creator_entry_number,
                                      'feature_name': site[0]
                                      }
                        Plasmid_Feature.add(self.tsession, input_dict)
                else:
                    input_dict = {'creator' : current_plasmid_entry.creator,
                                  'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                  'feature_name' : site[0]
                                  }
                    Plasmid_Feature.add(self.tsession, input_dict)


    def upload_file(self, current_plasmid_entry):
        buffy = self.generate_ape_file(current_plasmid_entry.get_id(), current_plasmid_entry.sequence, current_plasmid_entry.description, UID = current_plasmid_entry.plasmid_name)

        filename = current_plasmid_entry.plasmid_name or current_plasmid_entry.get_id()
        if isinstance(filename, str):
            filename = unicode(filename, 'utf-8')

        new_plasmid_file_entry = {'creator' : current_plasmid_entry.creator,
                                  'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                  'file_name' : filename,
                                  'file_type' : u'GenBank file',
                                  'Description' : current_plasmid_entry.description,
                                  'File' : buffy
                                  }

        Plasmid_File.add(self.tsession, new_plasmid_file_entry, silent=False)


    def generate_ape_file(self, Database_ID, complete_assembly, complete_description, UID = None,  mutations = None, mutant_feature_tuple = None):
        features_list = []
        features_query = self.tsession.query(Feature, Feature_Type).filter(Feature.Feature_type == Feature_Type.Feature_type)
        possible_features = [(record.Feature_name, record.Feature_type, record.Feature_sequence, color.color, record.description) for record, color in features_query]

        # Special cases for mutations
        if mutant_feature_tuple != None:
            possible_features.append(mutant_feature_tuple)

        for site in possible_features:
            target = site[2].upper()
            #Find forward sequences from features database
            if complete_assembly.find(target) != -1:
                # Treat type II restriction enzymes as special cases
                if site[1] == 'Restrxn Type II':
                    for instance in re.finditer(target, complete_assembly):
                        features_list.append(
                            SeqFeature(
                                CompoundLocation(
                                    [
                                        FeatureLocation(instance.start(), instance.end()),
                                        FeatureLocation(instance.start() + 7, instance.start() + 11)
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
                else:
                    for instance in re.finditer(target, complete_assembly):
                        features_list.append(
                            SeqFeature(
                                FeatureLocation(instance.start(), instance.end()),
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

            if complete_assembly.find(self.reverse_complement(target)) != -1:
                if site[1] == 'Restrxn Type II':
                    for instance in re.finditer(self.reverse_complement(target), complete_assembly):
                        features_list.append(
                            SeqFeature(
                                CompoundLocation(
                                    [
                                        FeatureLocation(instance.start(), instance.end()),
                                        FeatureLocation(instance.start() - 5, instance.start() - 1)
                                    ]
                                ),
                                type=site[1],
                                strand=-1,
                                qualifiers={"label": site[0],
                                            "ApEinfo_label": site[0],
                                            "ApEinfo_fwdcolor": site[3],
                                            "ApEinfo_revcolor": site[3],
                                            "ApEinfo_graphicformat": "arrow_data {{0 1 2 0 0 -1} {} 0}"
                                            }
                            )
                        )
                else:
                    for instance in re.finditer(self.reverse_complement(target), complete_assembly):
                        features_list.append(
                            SeqFeature(
                                FeatureLocation(instance.start(),instance.end()),
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

        # Check that UID is not greater than 15 characters long or null, otherwise defaults to Database ID
        if UID:
            if len(UID) > 15 or len(UID) == 0:
                plasmid_name = Database_ID
            else:
                plasmid_name = UID
        else:
            plasmid_name = Database_ID

        sequence = SeqRecord( Seq(complete_assembly,
                                  IUPAC.unambiguous_dna),
                              id=Database_ID, # ACCESSION AND VERSION
                              name=plasmid_name, # LOCUS
                              description=complete_description,
                              features = features_list
                              )

        buffy = StringIO.StringIO()
        SeqIO.write(sequence, buffy, 'gb')

        return buffy.getvalue()


    def mutant_check(self, temp_sequence):
        sequence_pass = True
        if any(rxn_site in temp_sequence for rxn_site in ['GGTCTC', 'GAGACC', 'CGTCTC', 'GAGACG']):
            sequence_pass = False
        return sequence_pass


    def generate_mutant_sequence(self, Plasmid_Feature_ID, mutation_list, write_to_file = False):
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
                       'V' : ['GTG', 'GTT','GTC', 'GTA'],
                       'W' : ['TGG'],
                       'Y' : ['TAT', 'TAC']
                       }

        feature_query = self.tsession.query(Plasmid_Feature, Feature).filter(Plasmid_Feature.ID == Plasmid_Feature_ID).filter(Plasmid_Feature.feature_name == Feature.Feature_name )

        for plasmid_feature, feature in feature_query:
            Actual_WT_feature_sequence = feature.Feature_sequence
            WT_feature_sequence = feature.Feature_sequence
            # print plasmid_feature.ID
            # print plasmid_feature.creator
            # print plasmid_feature.creator_entry_number
            # print feature.Feature_name
            # print feature.Feature_type
            # print feature.Feature_sequence

        CDS_mutant_constituents = []
        for point_mutant in mutation_list:

            if WT_feature_sequence[3 * point_mutant[1] - 3 : -(len(WT_feature_sequence) - 3 * point_mutant[1])] not in codon_table[point_mutant[0]]:
                raise Plasmid_Exception('User input target codon and sequence codon identities do not match!')

            mutant_codon_choices = codon_table[point_mutant[2]]

            for codon in mutant_codon_choices:
                temp_sequence = WT_feature_sequence[ : 3 * point_mutant[1] - 3 ] + codon + WT_feature_sequence[ - ( len(WT_feature_sequence) - ( 3 * point_mutant[1] ) ) : ]
                if self.mutant_check(temp_sequence):
                    CDS_mutant_constituents.append([point_mutant, codon])
                    break
                else:
                    continue
            if self.mutant_check(temp_sequence) != True:
                raise Plasmid_Exception("No suitable codons were found for %s%s%s! You're going to have to do this one by hand :P" %(point_mutant[0], point_mutant[1], point_mutant[2]))

            WT_feature_sequence = temp_sequence.upper()

        sequence_query = self.tsession.query(Plasmid_Feature, Feature, Feature_Type, Plasmid)\
            .filter(Plasmid_Feature.ID == Plasmid_Feature_ID)\
            .filter(Plasmid.creator == Plasmid_Feature.creator)\
            .filter(Plasmid.creator_entry_number == Plasmid_Feature.creator_entry_number)\
            .filter(Plasmid_Feature.feature_name == Feature.Feature_name)\
            .filter(Feature.Feature_type == Feature_Type.Feature_type)

        for plasmid_feature, feature, feature_type, plasmid in sequence_query:
            database_ID = plasmid.get_id()
            user_ID = plasmid.plasmid_name
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

        # Generate plasmid sequence with mutant CDS
        CDS_index_start = WT_plasmid_sequence.find(Actual_WT_feature_sequence)
        Mutant_plasmid_sequence = WT_plasmid_sequence[ : CDS_index_start ] + WT_feature_sequence + WT_plasmid_sequence[ -( len(WT_plasmid_sequence) - (CDS_index_start + len(WT_feature_sequence)) ) : ]

        # Annotate mutant codons
        Mutant_codon_list = []

        for mutation in mutation_list:

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

        mutant_genbank_file = self.generate_ape_file(database_ID,
                                                     Mutant_plasmid_sequence,
                                                     WT_feature_description,
                                                     UID = user_ID,
                                                     mutations=Mutant_codon_list,
                                                     mutant_feature_tuple=mutant_feature_tuple)

        with open('{0}_Mut.gb'.format(database_ID), 'w+b') as file:
            file.write(mutant_genbank_file)

        return Mutant_plasmid_sequence, CDS_mutant_constituents, mutant_genbank_file

    def add_mutant_to_db(self, Plasmid_Feature_ID, CDS_mutant_constituents, auto_commit=False):

        description_temp = []

        print CDS_mutant_constituents
        for mutation in CDS_mutant_constituents:
            print mutation
            description_temp.append(mutation[0][0] + str(mutation[0][1]) + mutation[0][2])

        new_CDS_Mutant_entry = {'Plasmid_Feature_ID' : Plasmid_Feature_ID,
                                'creator' : self.user_ID,
                                'Description' : '.'.join(description_temp)
                                }

        current_mutant_entry = CDS_Mutant.add(self.tsession, new_CDS_Mutant_entry, silent=False)

        for constituent in CDS_mutant_constituents:
            new_CDS_mutant_constituent = {'ID' : current_mutant_entry.ID,
                                          'mutation' : constituent[1],
                                          'position' : constituent[0][1],
                                          'wt_AA' : constituent[0][0],
                                          'mut_AA' : constituent[0][2]
                                          }

            CDS_Mutant_Constituent.add(self.tsession, new_CDS_mutant_constituent)

        if auto_commit:
            self.tsession.commit()

    def fetch_mutant_from_db(self, user_query):
        mutant_info = self.tsession.query(CDS_Mutant, CDS_Mutant_Constituent, Plasmid_Feature, Plasmid).filter(
            and_(
                CDS_Mutant.ID == user_query,
                CDS_Mutant.ID == CDS_Mutant_Constituent.ID,
                Plasmid_Feature.ID == CDS_Mutant.Plasmid_Feature_ID,
                Plasmid_Feature.creator == Plasmid.creator,
                Plasmid_Feature.creator_entry_number == Plasmid.creator_entry_number
            )
        )

        # mutant_info = self.tsession.execute('''
        #     SELECT CDS_Mutant.ID AS CDS_Mutant_ID, CDS_Mutant.Plasmid_Feature_ID AS CDS_Mutant_Plasmid_Feature_ID, CDS_Mutant.Mutant_ID AS CDS_Mutant_Mutant_ID, CDS_Mutant.creator AS CDS_Mutant_creator, CDS_Mutant.Description AS CDS_Mutant_Description, CDS_Mutant.date AS CDS_Mutant_date, CDS_Mutant_Constituent.ID AS CDS_Mutant_Constituent_ID, CDS_Mutant_Constituent.mutation AS CDS_Mutant_Constituent_mutation, CDS_Mutant_Constituent.position AS CDS_Mutant_Constituent_position, CDS_Mutant_Constituent.wt_AA AS CDS_Mutant_Constituent_wt_AA, CDS_Mutant_Constituent.mut_AA AS CDS_Mutant_Constituent_mut_AA,
        #     CDS_Mutant_Constituent.Description AS CDS_Mutant_Constituent_Description, Plasmid_Feature.ID AS Plasmid_Feature_ID, Plasmid_Feature.creator AS Plasmid_Feature_creator, Plasmid_Feature.creator_entry_number AS Plasmid_Feature_creator_entry_number,
        #     Plasmid_Feature.feature_name AS Plasmid_Feature_feature_name, Plasmid.creator AS Plasmid_creator, Plasmid.creator_entry_number AS Plasmid_creator_entry_number, Plasmid.plasmid_name AS Plasmid_plasmid_name, Plasmid.plasmid_type AS Plasmid_plasmid_type, Plasmid.location AS Plasmid_location,
        #     Plasmid.description AS Plasmid_description, Plasmid.sequence AS Plasmid_sequence, Plasmid.status AS Plasmid_status, Plasmid.date AS Plasmid_date
        #     FROM CDS_Mutant
        #     INNER JOIN CDS_Mutant_Constituent ON CDS_Mutant.ID = CDS_Mutant_Constituent.ID
        #     INNER JOIN Plasmid_Feature ON Plasmid_Feature.ID = CDS_Mutant.Plasmid_Feature_ID
        #     INNER JOIN Plasmid ON Plasmid_Feature.creator = Plasmid.creator AND Plasmid_Feature.creator_entry_number = Plasmid.creator_entry_number
        #     WHERE CDS_Mutant.ID = :user_query
        #   ''', dict(user_query=user_query))

        # SELECT * FROM CDS_Mutant INNER JOIN CDS_Mutant_Constituent ON CDS_Mutant.ID == CDS_Mutant_Constituent.ID

        mutation_tuples = []
        for cds_mutant, cds_mutant_constituent, plasmid_feature, plasmid in mutant_info:
            print cds_mutant_constituent.position
            Plasmid_Feature_ID = cds_mutant.Plasmid_Feature_ID
            mutation_tuples.append((cds_mutant_constituent.wt_AA, cds_mutant_constituent.position, cds_mutant_constituent.mut_AA))
            database_ID = plasmid.get_id()
            mutant_ID = cds_mutant.Mutant_ID

        Mutant_plasmid_sequence, CDS_mutant_constituents, mutant_genbank_file = self.generate_mutant_sequence(Plasmid_Feature_ID, mutation_tuples, write_to_file=True)

        with open('{0}_MutID-{1:04d}.gb'.format(database_ID, mutant_ID), 'w+b') as file:
            file.write(mutant_genbank_file)

    def generate_ape_from_database_ID(self, creator, creator_entry_number, write_to_file = False):
        my_plasmid = self.tsession.query(Plasmid).filter(Plasmid.creator == creator).filter(Plasmid.creator_entry_number == creator_entry_number).one()
        database_ID = my_plasmid.get_id()
        genbank_file = self.generate_ape_file(database_ID, my_plasmid.sequence.upper(), my_plasmid.description, UID = my_plasmid.plasmid_name)

        if write_to_file == True:
            with open('{0}.gb'.format(database_ID), 'w+b') as file:
                file.write(genbank_file)

        return genbank_file

    def update_features(self):
        all_my_plasmids = self.tsession.query(Plasmid)
        for plasmid in all_my_plasmids:
            self.add_features(plasmid)
        self.tsession.commit()
        print 'Update Complete!'


    def nuke_plasmids(self, nuke_list):
        """
        Deletes plasmid entries from the database
        :param nuke_list: list of (creator, creator_entry_number) tuples of plasmids to be deleted
        :return:
        """
        for target in nuke_list:
            try:
                target_query = self.tsession.query(Plasmid).filter(and_(Plasmid.creator == target[0], Plasmid.creator_entry_number == target[1])).one()
            except:
                print "p{0}{1:04d} doesn't seem to exist. Moving on!".format(target[0], target[1])
                continue

            # Checking for special cases where a plasmid might be particularly important and we won't want to delete it
            mutant_query = self.tsession.query(Plasmid_Feature, CDS_Mutant).filter(and_(Plasmid_Feature.creator == target[0],
                                                                                             Plasmid_Feature.creator_entry_number == target[1],
                                                                                             Plasmid_Feature.ID == CDS_Mutant.Plasmid_Feature_ID)).count()
            publication_query = self.tsession.query(Publication_Plasmid).filter(and_(Publication_Plasmid.creator == target[0], Publication_Plasmid.creator_entry_number == target[1])).count()

            if mutant_query != 0:
                print "WARNING: A database record in CDS_Mutants references p{0}{1:04d}. This plasmid will not be deleted.".format(target[0], target[1])
                continue
            if publication_query != 0:
                print "WARNING: A database record in Publication_Plasmid references p{0}{1:04d}. This plasmid will not be deleted.".format(target[0], target[1])
                continue

            if target_query.plasmid_type == 'part':
                try:
                    dsession = self.dbi.get_session()
                    dsession.query(Part_Plasmid_Part).filter(and_(Part_Plasmid_Part.creator == target[0], Part_Plasmid_Part.creator_entry_number == target[1])).delete()
                    dsession.query(Part_Plasmid).filter(and_(Part_Plasmid.creator == target[0], Part_Plasmid.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_Feature).filter(and_(Plasmid_Feature.creator == target[0], Plasmid_Feature.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_File).filter(and_(Plasmid_File.creator == target[0], Plasmid_File.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid).filter(and_(Plasmid.creator == target[0], Plasmid.creator_entry_number == target[1])).delete()
                    dsession.commit()
                    dsession.close()
                    print '***p{0}{1:04d} DELETED***'.format(target[0], target[1])
                except Exception, e:
                    import traceback

                    print('Failure')
                    print(str(e))
                    print(traceback.format_exc())

                    if dsession:
                        dsession.rollback()
                        dsession.close()

            elif target_query.plasmid_type == 'cassette':
                try:
                    dsession = self.dbi.get_session()
                    dsession.query(Cassette_Assembly).filter(and_(Cassette_Assembly.Cassette_creator == target[0],Cassette_Assembly.Cassette_creator_entry_number == target[1])).delete()
                    dsession.query(Cassette_Plasmid).filter(and_(Cassette_Plasmid.creator == target[0],Cassette_Plasmid.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_Feature).filter(and_(Plasmid_Feature.creator == target[0],Plasmid_Feature.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_File).filter(and_(Plasmid_File.creator == target[0],Plasmid_File.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid).filter(and_(Plasmid.creator == target[0], Plasmid.creator_entry_number == target[1])).delete()
                    dsession.commit()
                    dsession.close()
                    print '***p{0}{1:04d} DELETED***'.format(target[0], target[1])
                except Exception, e:
                    import traceback

                    print('Failure')
                    print(str(e))
                    print(traceback.format_exc())

                    if dsession:
                        dsession.rollback()
                        dsession.close()

            elif target_query.plasmid_type == 'multicassette':
                try:
                    dsession = self.dbi.get_session()
                    dsession.query(Multicassette_Assembly).filter(and_(Multicassette_Assembly.Multicassette_creator == target[0], Multicassette_Assembly.Multicassette_creator_entry_number == target[1])).delete()
                    dsession.query(Multicassette_Plasmid).filter(and_(Multicassette_Plasmid.creator == target[0], Multicassette_Plasmid.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_Feature).filter(and_(Plasmid_Feature.creator == target[0], Plasmid_Feature.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_File).filter(and_(Plasmid_File.creator == target[0], Plasmid_File.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid).filter(and_(Plasmid.creator == target[0], Plasmid.creator_entry_number == target[1])).delete()
                    dsession.commit()
                    dsession.close()
                    print '***p{0}{1:04d} DELETED***'.format(target[0], target[1])
                except Exception, e:
                    import traceback

                    print('Failure')
                    print(str(e))
                    print(traceback.format_exc())

                    if dsession:
                        dsession.rollback()
                        dsession.close()

            elif target_query.plasmid_type == 'other':
                try:
                    dsession = self.dbi.get_session()
                    dsession.query(Plasmid_Feature).filter(and_(Plasmid_Feature.creator == target[0],Plasmid_Feature.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_File).filter(and_(Plasmid_File.creator == target[0],Plasmid_File.creator_entry_number == target[1])).delete()
                    dsession.query(Other_Plasmid).filter(and_(Other_Plasmid.creator == target[0], Other_Plasmid.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid).filter(and_(Plasmid.creator == target[0], Plasmid.creator_entry_number == target[1])).delete()
                    dsession.commit()
                    dsession.close()
                    print '***p{0}{1:04d} DELETED***'.format(target[0], target[1])
                except Exception, e:
                    import traceback

                    print('Failure')
                    print(str(e))
                    print(traceback.format_exc())

                    if dsession:
                        dsession.rollback()
                        dsession.close()

            elif target_query.plasmid_type == 'template':
                try:
                    dsession = self.dbi.get_session()
                    dsession.query(Plasmid_Feature_Design).filter(and_(Plasmid_Feature_Design.parent_creator == target[0],Plasmid_Feature_Design.parent_creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_Feature).filter(and_(Plasmid_Feature.creator == target[0],Plasmid_Feature.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_File).filter(and_(Plasmid_File.creator == target[0],Plasmid_File.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid).filter(and_(Plasmid.creator == target[0], Plasmid.creator_entry_number == target[1])).delete()
                    dsession.commit()
                    dsession.close()
                    print '***p{0}{1:04d} DELETED***'.format(target[0], target[1])
                except Exception, e:
                    import traceback

                    print('Failure')
                    print(str(e))
                    print(traceback.format_exc())

                    if dsession:
                        dsession.rollback()
                        dsession.close()

            elif target_query.plasmid_type == 'design':
                try:
                    dsession = self.dbi.get_session()
                    dsession.query(Plasmid_Feature_Design).filter(and_(Plasmid_Feature_Design.child_creator == target[0],Plasmid_Feature_Design.child_creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_Feature).filter(and_(Plasmid_Feature.creator == target[0],Plasmid_Feature.creator_entry_number == target[1])).delete()
                    dsession.query(Design_Plasmid).filter(and_(Design_Plasmid.creator == target[0],Design_Plasmid.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid_File).filter(and_(Plasmid_File.creator == target[0],Plasmid_File.creator_entry_number == target[1])).delete()
                    dsession.query(Plasmid).filter(and_(Plasmid.creator == target[0], Plasmid.creator_entry_number == target[1])).delete()
                    dsession.commit()
                    dsession.close()
                    print '***p{0}{1:04d} DELETED***'.format(target[0], target[1])
                except Exception, e:
                    import traceback

                    print('Failure')
                    print(str(e))
                    print(traceback.format_exc())

                    if dsession:
                        dsession.rollback()
                        dsession.close()

            else:
                raise Plasmid_Exception("Plasmid type not recognized... What kind of plasmid is this?!?!?!?!")


    def design_feature(self, design_list, designed_plasmid_ID, designed_plasmid_description, designed_plasmid_location, tsession = None):
        """
        This function will allow users to replace features in annotated template plasmids with their own designs
        :param design_list: a list of dicts to produce a new plasmid where each dict has the following keys/values:
                feature_ID: Feature ID found in Plasmid_Feature.ID
                feature_design_sequence: Sequence of the designed feature
                feature_design_name: Name of the designed feature
                feature_design_description: Description of the designed feature
        :param designed_plasmid_ID: Essentially the User ID for this designed plasmid
        :param designed_plasmid_description: Description for the designed plasmid
        :param designed_plasmid_location: Design location...
        :return:
        """

        # The user should handle transaction controls on any session passed in
        auto_commit = False
        if not tsession:
            auto_commit = True
            tsession = self.tsession

        WT_set = tsession.query(Plasmid_Feature, Feature, Plasmid).filter(
            and_(Plasmid_Feature.feature_name == Feature.Feature_name,
                 Plasmid.creator == Plasmid_Feature.creator,
                 Plasmid.creator_entry_number == Plasmid_Feature.creator_entry_number))

        # Make sure all features are from the same plasmid
        plasmid_query_list = [WT_set.filter(Plasmid_Feature.ID == design['feature_ID']).one() for design in design_list]

        parent_plasmid_set = set([(feature_source.Plasmid.creator, feature_source.Plasmid.creator_entry_number) for feature_source in plasmid_query_list])

        if len(parent_plasmid_set) != 1:
            raise Plasmid_Exception('All target features need to be from the same Plasmid!')

        ####################
        # Generate designs #
        ####################

        # Get information on parent plasmid
        parent_creator = plasmid_query_list[0].Plasmid.creator
        parent_creator_entry_number = plasmid_query_list[0].Plasmid.creator_entry_number
        WT_sequence = plasmid_query_list[0].Plasmid.sequence.upper()
        WT_plasmid_type = plasmid_query_list[0].Plasmid.plasmid_type
        design_intermediate = WT_sequence.upper()

        for feature_design in design_list:
            feature_query = WT_set.filter(Plasmid_Feature.ID == feature_design['feature_ID']).one()
            WT_feature = feature_query.Feature.Feature_sequence

            # Adding support for features that occur more than once in a plasmid
            # If instances > 1, then replace features in reverse so that feature indicies from WT remain the same

            if len([instance for instance in re.finditer(WT_feature.upper(), design_intermediate.upper())]) == 1:
                for instance in re.finditer(WT_feature.upper(), design_intermediate.upper()):
                    design_intermediate = design_intermediate[:instance.start()] + feature_design['feature_design_sequence'].upper() + design_intermediate[instance.end():]

            elif len([instance for instance in re.finditer(WT_feature.upper(), design_intermediate.upper())]) > 1:
                feature_indicies = [{'start': instance.start(), 'end': instance.end()} for instance in re.finditer(WT_feature.upper(), design_intermediate.upper())]

                for feature_index in reversed(feature_indicies):
                    design_intermediate = design_intermediate[:feature_index['start']] + feature_design['feature_design_sequence'].upper() + design_intermediate[feature_index['end']:]

            else:
                raise Plasmid_Exception('An instance of {0} was not found in p{1}{2:04d}!'.format(feature_query.Plasmid_Feature.feature_name,
                                                                                                  feature_query.Plasmid.creator,
                                                                                                  feature_query.Plasmid.creator_entry_number))

            # todo: push to Features if the new designed feature does not exist, otherwise notify the user and use the already existing feature
            # Add designed feature to Features if it does not exist. Otherwise notify the user and use the already existing feature.

            all_features_query = tsession.query(Feature)
            existing_features_check = set(feature.Feature_sequence.upper() for feature in all_features_query)
            if feature_design['feature_design_sequence'].upper() in existing_features_check:
                feature_already_exists = all_features_query.filter(Feature.Feature_sequence == feature_design['feature_design_sequence']).one()
                print "Your designed feature \"{0}\" (F) already exists in the Features database as \"{1}\"".format(feature_design['feature_design_name'], feature_already_exists.Feature_name)
            elif self.reverse_complement(feature_design['feature_design_sequence'].upper()) in existing_features_check:
                feature_already_exists = all_features_query.filter(Feature.Feature_sequence == feature_design['feature_design_sequence']).one()
                print "Your designed feature \"{0}\" (R) already exists in the Features database as \"{1}\"".format(feature_design['feature_design_name'], feature_already_exists.Feature_name)
            else:
                feature_input_dict = {'Feature_name': feature_design['feature_design_name'],
                                      'Feature_type': feature_query.Feature.Feature_type,
                                      'Feature_sequence': feature_design['feature_design_sequence'],
                                      'description': feature_design['feature_design_sequence']}
                Feature.add(tsession, feature_input_dict)

        final_designed_sequence = design_intermediate

        # Implement plasmid_checks()
        input_dict = {'sequence': final_designed_sequence,
                      'plasmid_name': designed_plasmid_ID or ' with designed feature \"{0}\" '.format(feature_design['feature_design_name'])
                      }

        if WT_plasmid_type == 'part':
            plasmid_parts_query = self.tsession.query(Part_Plasmid_Part).filter(and_(
                Part_Plasmid_Part.creator == parent_creator,
                Part_Plasmid_Part.creator_entry_number == parent_creator_entry_number
            ))
            plasmid_parts_string = ', '.join([part_part.part_number for part_part in plasmid_parts_query])

            self.plasmid_checks(input_dict, WT_plasmid_type, plasmid_parts_string)
            print plasmid_parts_string
        else:
            self.plasmid_checks(input_dict, WT_plasmid_type)

        ###############################
        # Push things to the database #
        ###############################

        plasmid_input_dict = {'creator': self.user_ID,
                              'plasmid_name': designed_plasmid_ID,
                              'plasmid_type': WT_plasmid_type,
                              'location': designed_plasmid_location,
                              'description': designed_plasmid_description,
                              'sequence': final_designed_sequence,
                              'status': u'designed'
                              }

        #Push Design_Plasmid
        design_Plasmid_entry = Plasmid.add(tsession, plasmid_input_dict, silent=False)

        # Check if the parent plasmid was the result of a previous round of design
        # If so, use the same common ancestor; else, set the common ancestor to the parent plasmid
        plasmid_feature_design_query = self.tsession.query(Plasmid_Design_Ancestors).filter(
            Plasmid_Design_Ancestors.child_creator == list(parent_plasmid_set)[0][0],
            Plasmid_Design_Ancestors.child_creator_entry_number == list(parent_plasmid_set)[0][1]
        ).first()

        if plasmid_feature_design_query is None:
            plasmid_design_ancestors_input_dict = {'parent_creator': list(parent_plasmid_set)[0][0],
                                                   'parent_creator_entry_number': list(parent_plasmid_set)[0][1],
                                                   'child_creator': design_Plasmid_entry.creator,
                                                   'child_creator_entry_number': design_Plasmid_entry.creator_entry_number,
                                                   'common_ancestor_creator': list(parent_plasmid_set)[0][0],
                                                   'common_ancestor_creator_entry_number': list(parent_plasmid_set)[0][
                                                       1]
                                                   }
            Plasmid_Design_Ancestors.add(tsession, plasmid_design_ancestors_input_dict, silent=False)
        else:
            plasmid_design_ancestors_input_dict = {'parent_creator': list(parent_plasmid_set)[0][0],
                                                   'parent_creator_entry_number': list(parent_plasmid_set)[0][1],
                                                   'child_creator': design_Plasmid_entry.creator,
                                                   'child_creator_entry_number': design_Plasmid_entry.creator_entry_number,
                                                   'common_ancestor_creator': plasmid_feature_design_query.common_ancestor_creator,
                                                   'common_ancestor_creator_entry_number': plasmid_feature_design_query.common_ancestor_creator_entry_number,
                                                   }
            Plasmid_Design_Ancestors.add(tsession, plasmid_design_ancestors_input_dict, silent=False)


        for feature_design in design_list:
            feature_design_input_dict = {'design_creator': design_Plasmid_entry.creator,
                                         'design_creator_entry_number': design_Plasmid_entry.creator_entry_number,
                                         'feature_ID': feature_design['feature_ID']
                                         }

            feature_design_entry = Plasmid_Design_Features.add(tsession, feature_design_input_dict, silent=False)

        self.add_features(design_Plasmid_entry)
        
        if auto_commit:
            tsession.commit()

        return design_Plasmid_entry


    def return_designable_features(self, engine, creator, creator_entry_number):
        '''
        Returns a list of designable features for a given plasmid.

        :param engine: An SQLAlchemy engine object e.g. engine = dbi.get_engine()
        :param creator: The UserID of the creator e.g. "JL"
        :param creator_entry_number: The entry number of the plasmid in the user's namespace e.g. 1
        :return: A list of dicts containing attributes of all designable features for the plasmid. Attributes include:  feature_ID, creator, creator_entry_number, feature_name, MD5_hash, sequence, description, feature_type, color.
        '''

        return call_procedure(engine, 'getDesignableFeatures', [creator, creator_entry_number], as_dict = True)

    def return_plasmid_resistance(self, creator, creator_entry_number):
        """
        Returns all resistances for a given plasmid as derived from the Plasmid_Features table

        :param creator: creator
        :param creator_entry_number: creator_entry_number
        :return:
        """
        plasmid_resistances = self.tsession.execute("SELECT * FROM Feature INNER JOIN Plasmid_Feature on Feature.Feature_name = Plasmid_Feature.feature_name where Feature.Feature_type = 'Resistance' and Plasmid_Feature.creator = :creator_ and Plasmid_Feature.creator_entry_number = :creator_entry_number_",
        {"creator_" : creator, "creator_entry_number_" : creator_entry_number})

        plasmid_resistance_list = []

        # All resistance features are of the form {RESISTANCE (Source)} e.g. {CARB (Kale)}
        # Doing .split()[0] returns just the resistance name as it appears in the resistance table e.g. CARB
        for plasmid_resistance in plasmid_resistances:
            plasmid_resistance_list.append(dict(plasmid_resistance)['feature_name'].split()[0])

        # List comprehension doesn't want to work for me for some reason...
        # plasmid_resistance_list = [dict(plasmid_resistance)['feature_name'].split()[0] for plasmid_resistance in plasmid_resistances]

        print sorted(plasmid_resistance_list)
        return sorted(plasmid_resistance_list)

    def add_new_user_defined_features(self, feature_dict, auto_commit=False):
        """
        This function will allow users to add new featurse to the Feature table given that they satisfy the following criteria:
        1. New features cannot be a subset of an existing feature (MAKE SURE TO CHECK F AND R!)
        2. New features cannot contain an existing feature
        3. New features should be at least ~15bp long. This is an arbitrary cutoff, but shorter sequences have a greater
           chance of popping up randomly and causing problems later down the road
        4. Feature Names must be unique
        5. I'll think of more stuff as this is written...
        :param feature_dict: dict with keys "feature_name", "feature_type", "feature_description", and "feature_sequence"
        :return: New Feature SQLAlchemy object
        """
        import re

        # Get Feature information
        database_features = self.tsession.query(Feature)

        # Create lists of database features and sequences
        database_feature_names = [database_feature.Feature_name.lower() for database_feature in database_features]
        database_feature_sequences = [database_feature.Feature_sequence.upper() for database_feature in database_features]
        # Remove restriction sites from sequence list
        database_feature_sequences.remove('CGTCTC')
        database_feature_sequences.remove('GGTCTC')

        # Test print
        import pprint
        pprint.pprint(database_feature_names)
        pprint.pprint(database_feature_sequences)


        ##################################
        # Perform Checks on New Features #
        ##################################
        # todo: add naming check for resistance features in the form of "[Standard Resistance Abbreviation] [Descriptor]"
        # New features should be at least ~15bp long
        if len(feature_dict['Feature_sequence']) < 15:
            raise Plasmid_Exception('Features must be 15bp or longer!')
        for database_feature in database_features:
            # Feature names must be unique
            if database_feature.Feature_name.lower() == feature_dict['Feature_name']:
                raise Plasmid_Exception('Feature names must be unique!')
            # New features cannot contain other features (F)
            if re.search(database_feature.Feature_sequence.upper(), feature_dict['Feature_sequence'].upper()) != None:
                raise Plasmid_Exception('{0} contains a subset sequence from {1}'.format(feature_dict['Feature_name'], database_feature.Feature_name))
            # New features cannot contain other features (R)
            if re.search(self.reverse_complement(database_feature.Feature_sequence.upper()), feature_dict['Feature_sequence'].upper()) != None:
                raise Plasmid_Exception('{0} contains a subset sequence from {1}'.format(feature_dict['Feature_name'], database_feature.Feature_name))
            # Existing features cannot contain a new feature (F)
            if re.search(feature_dict['Feature_sequence'].upper(), database_feature.Feature_sequence.upper()) != None:
                raise Plasmid_Exception('An existing feature ({0}) contains the sequence for {1}'.format(database_feature.Feature_name, feature_dict['Feature_name']))
            # Existing features cannot contain a new feature (R)
            if re.search(self.reverse_complement(feature_dict['Feature_sequence'].upper()), database_feature.Feature_sequence.upper()) != None:
                raise Plasmid_Exception('An existing feature ({0}) contains the sequence for {1}'.format(database_feature.Feature_name, feature_dict['Feature_name']))

        # Commit to Database
        new_feature = Feature.add(self.tsession, feature_dict, silent=True)

        if auto_commit:
            self.tsession.commit()

        return new_feature