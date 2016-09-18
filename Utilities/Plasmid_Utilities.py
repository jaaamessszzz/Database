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
    from db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature, Plasmid_Feature, Plasmid_File, CDS_Mutant, CDS_Mutant_Constituent
except:
    # nasty hack since we are not packaging things up properly yet for external use (e.g. the website)
    from kprimers.db.interface import DatabaseInterface
    from kprimers.db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature, Plasmid_Feature, Plasmid_File, CDS_Mutant, CDS_Mutant_Constituent

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

    ###########################
    # Planned Funtionality:   #
    ###########################

    * Multicassette assembly and submission
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
            self.username = os.getlogin() # The webserver cannot call this function so it must instead pass username in as an argument
        self.user_ID = self.tsession.query(Users).filter(Users.lab_username == self.username).one().ID


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

        # Start with Part one (cassette assembly) so that there is consistency in how the plasmids are assembled
        if assembly_type.lower() == 'cassette':
            table_info = self.fetch_cassette_parts(input_sequences, table_info)
            for part in table_info['Part list']:
                if part[:4].upper() == 'CCCT':
                    intermediate = part.upper()
                    table_info['Part list'].remove(part)
                    break

        if assembly_type.lower() == 'part':
            # Checks that part 3s have len % 3 = 0 and do not contain a starting MET or stop codon anywhere
            if '3' in part_type:
                if len(input_sequences[0]) % 3 != 0:
                    raise Plasmid_Exception('Your input sequence contains incomplete codons!!')

            if any(input_type == part_type for input_type in ['3', '3a']):
                if input_sequences[0][:3].upper() == 'ATG':
                    # Do we need to raise? Or can we just remove the start codon for the user and inform them?
                    raise Plasmid_Exception(
                        'Please remove the start codon from your input sequence! The start codon is already in the part plasmid sequence.')

                # Check that stop codons are not in the input coding sequence
                if any(codon in [input_sequences[0][i:i + 3].upper() for i in range(0, len(input_sequences[0]), 3)] for codon in ['TAA', 'TAG', 'TGA']):
                    raise Plasmid_Exception('There are stop codons in your coding sequence! Plase remove them!')

            # Stop codons: TAA TAG TGA

            part_entry_vector = self.tsession.query(Plasmid).filter(Plasmid.creator == 'JL').filter(Plasmid.creator_entry_number == 2).one()
            sequence_upper = self.add_part_arms(input_sequences[0], part_type)
            table_info['Part list'].append(sequence_upper[ sequence_upper.find('CGTCTC') + 7 : sequence_upper.find('GAGACG') - 1])
            intermediate = part_entry_vector.sequence[ : (part_entry_vector.sequence.upper().find('GAGACG') - 1) ]

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

        if intermediate[-4:] != intermediate[:4]:
            raise Plasmid_Exception('Incomplete assembly! This assembly does not produce a circular plasmid! :(')

        table_info['Complete Assembly'] = intermediate[:-4].upper()
        table_info['Complete Description'] = ' | '.join(table_info['Description'])

        return table_info

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
                raise Plasmid_Exception('There is more than one forward %s site in %s! %s' % (restriction_enzyme, plasmid.plasmid_name, plasmid.sequence.count(site_F)))
            if sequence_upper.count(site_R) != 1:
                raise Plasmid_Exception('There is more than one reverse %s site in %s!' % (restriction_enzyme, plasmid.plasmid_name))

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
            part_sequences.append(sequence_upper[sequence_upper.find(left_overhang, sequence_upper.find(site_F)):sequence_upper.find(right_overhang,sequence_upper.find(site_R) - 10) + 4])

        if table_info:
            return table_info
        else:
            return part_sequences

    def add_part_arms(self, sequence, part_type):
        # Some of the overhangs aren't added yet since I don't have access to APE at the moment...
        arm_dict = {'1'  : ['', ''],
                    '2'  : ['gcatCGTCTCaAGCAGGTCTCAAACG', 'TATGtGAGACCtGAGACGgcat'],
                    # '2a' : ['', ''],
                    #  '2b' : ['', ''],
                    '3'  : ['gcatCGTCTCaAGCAGGTCTCaTATG', 'ATCCtGAGACCtGAGACGgcat'],
                    '3a' : ['gcatCGTCTCaAGCAGGTCTCaTATG', 'GGTAGCGGCAGCGGCAGCTTCTtGAGACCtGAGACGgcat'],
                    '3b' : ['gcatCGTCTCaAGCAGGTCTCATTCT', 'ATCCtGAGACCtGAGACGgcat'],
                    '4'  : ['gcatCGTCTCaAGCAGGTCTCaATCC', 'GCTGtGAGACCtGAGACGgcat'],
                    # '4a' : ['gcatCGTCTCaAGCAGGTCTCaATCC', ''],
                    #  '4b' : ['', 'GCTGtGAGACCtGAGACGgcat'],
                    '5'  : ['', ''],
                    '6'  : ['gcatCGTCTCaAGCAGGTCTCA', 'TGAGACCtGAGACGgcat'],
                    '7'  : ['', '']
                    }

        try:
            part_sequence = arm_dict[part_type][0] + sequence + arm_dict[part_type][1]
            if part_sequence.upper().find('CGTCTC') == -1 or part_sequence.upper().find('GAGACG') == -1:
                raise Plasmid_Exception('BsmBI sites not found! Something fishy is going on...')

        except:
            raise Plasmid_Exception('%s part types are not supported yet!' %part_type)

        return part_sequence.upper()


    def plasmid_checks(self, input_dict, assembly_type, part_type = None):
        '''
        Verifies that a plasmid conforms to modular cloning sequence standards before entry into database

        # Part Plasmids #
          * Part 1 and Part 5 should contain one reverse BsmBI site
          * If not a Part 1 or Part 5, sequence should not contain BsmBI sites
          *
        :return:
        '''
        if assembly_type.lower() == 'part':
            if part_type == '1' or part_type == '5':
                if input_dict['sequence'].count('GAGACG') != 1:
                    raise Plasmid_Exception('There should be only one reverse %s site in a Connector Part! Your part %s submission contains %s reverse %s sites.' % (
                    'BsmBI', part_type, input_dict['sequence'].count('GAGACG'), 'BsmBI'))
            else:
                if input_dict['sequence'].count('GAGACG') != 0:
                    raise Plasmid_Exception('There should not be any reverse %s sites in this part type! Your submission %s contains %s reverse %s site(s).' % (
                        'BsmBI', input_dict['plasmid_name'], input_dict['sequence'].count('GAGACG'), 'BsmBI'))
            if input_dict['sequence'].count('CGTCTC') != 0:
                raise Plasmid_Exception('There is more than one forward %s site in %s! This plasmid contains %s forward %s sites.' % (
                    'BsmBI', input_dict['plasmid_name'], input_dict['sequence'].count('CGTCTC'), 'BsmBI'))
            if input_dict['sequence'].count('GGTCTC') != 1:
                raise Plasmid_Exception('There is more than one forward %s site in %s! Your submission contains %s forward %s sites.' % (
                    'BsaI', input_dict['plasmid_name'], input_dict['sequence'].count('GGTCTC'), 'BsaI'))
            if input_dict['sequence'].count('GAGACC') != 1:
                raise Plasmid_Exception('There is more than one reverse %s site in %s! Your submission contains %s reverse %s sites.' % (
                    'BsaI', input_dict['plasmid_name'], input_dict['sequence'].count('GAGACC'), 'BsaI'))

        if assembly_type.lower() == 'cassette':
            if input_dict['sequence'].count('GGTCTC') != 0:
                raise Plasmid_Exception('There are BsaI sites in your completed cassette assembly! This thing is going to get wrecked during the final digestion step of the Golden Gate Assembly!')
            if input_dict['sequence'].count('GAGACC') != 0:
                raise Plasmid_Exception('There are BsaI sites in your completed cassette assembly! This thing is going to get wrecked during the final digestion step of the Golden Gate Assembly!')

            # I need to implement more complicated logic when checking for cassete connector sites in cassette assemblies
            # It looks like LS contains a forward BsmBI site and all other connectors have a reverse BsmBI site...
            # Not going to deal with that right now...

            # if input_dict['sequence'].count('CGTCTC') != 0:
            #     raise Plasmid_Exception('Your submission contains %s forward %s sites. There should be none!' % ( input_dict['sequence'].count('CGTCTC'), 'BsmBI'))
            # if input_dict['sequence'].count('GAGACG') != 2:
            #     raise Plasmid_Exception('There should only be two reverse %s sites in a cassette assembly. Your submission contains %s reverse %s sites.' % ('BsmBI', input_dict['sequence'].count('GAGACG'), 'BsmBI'))

        if assembly_type.lower() == 'multicassette':
            # Ehhh... I'll leave this for now. There really aren't any restrictions for a multicasssette plasmid that I can think of for now
            if input_dict['sequence'].count('GGTCTC') != 1:
                raise Plasmid_Exception('There is more than one forward %s site in %s! Your submission contains %s forward %s sites.' % (
                'BsaI', input_dict['plasmid_name'], input_dict['sequence'].count('GGTCTC'), 'BsaI'))
            if input_dict['sequence'].count('GAGACC') != 1:
                raise Plasmid_Exception('There is more than one reverse %s site in %s! Your submission contains %s reverse %s sites.' % (
                'BsaI', input_dict['plasmid_name'], input_dict['sequence'].count('GAGACC'), 'BsaI'))


    def add_part_plasmid_to_db(self, input_dict, part_type, auto_commit = True):
        '''
        :param user_input:
        :param table_info:
        :param auto_commit:
        :return: The created Plasmid object.
        '''

        current_plasmid_entry = Plasmid.add(self.tsession, input_dict, silent=True)

        new_part_plasmid_entry = {'creator' : current_plasmid_entry.creator,
                                  'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                  'resistance' : 'CM'
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


    def add_cassette_plasmid_to_db(self, input_dict, table_info, auto_commit = True):

        current_plasmid_entry = Plasmid.add(self.tsession, input_dict, silent=False)

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


    def add_features(self, current_plasmid_entry):
        features_query = self.tsession.query(Feature, Feature_Type).filter(
            Feature.Feature_type == Feature_Type.Feature_type)
        possible_features = [
            (record.Feature_name, record.Feature_type, record.Feature_sequence, color.color, record.description) for
            record, color in features_query]
        for site in possible_features:
            target = site[2].upper()
            if current_plasmid_entry.sequence.find(target) != -1 or current_plasmid_entry.sequence.find(self.reverse_complement(target)) != -1:
                input_dict = {'creator' : current_plasmid_entry.creator,
                              'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                              'feature_name' : site[0]
                              }
                Plasmid_Feature.add(self.tsession, input_dict)


    def upload_file(self, current_plasmid_entry):
        buffy = self.generate_ape_file(current_plasmid_entry.get_id(), current_plasmid_entry.sequence, current_plasmid_entry.description, UID = current_plasmid_entry.plasmid_name)

        filename = current_plasmid_entry.plasmid_name or current_plasmid_entry.get_id()
        new_plasmid_file_entry = {'creator' : current_plasmid_entry.creator,
                                  'creator_entry_number' : current_plasmid_entry.creator_entry_number,
                                  'file_name' : unicode(filename, 'utf-8'),
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
                       'V' : ['GTG', 'GTT' ,'GTC', 'GTA'],
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
        mutant_info = self.tsession.query(CDS_Mutant, CDS_Mutant_Constituent, Plasmid_Feature, Plasmid)\
            .filter(CDS_Mutant.ID == CDS_Mutant_Constituent.ID)\
            .filter(Plasmid_Feature.ID == CDS_Mutant.Plasmid_Feature_ID)\
            .filter(Plasmid_Feature.creator == Plasmid.creator)\
            .filter(Plasmid_Feature.creator_entry_number == Plasmid.creator_entry_number)\
            .filter(CDS_Mutant.ID == user_query)

        # SELECT * FROM CDS_Mutant INNER JOIN CDS_Mutant_Constituent ON CDS_Mutant.ID == CDS_Mutant_Constituent.ID

        mutation_tuples = []
        for cds_mutant, cds_mutant_constituent, plasmid_feature, plasmid in mutant_info:
            Plasmid_Feature_ID = cds_mutant.Plasmid_Feature_ID
            mutation_tuples.append((cds_mutant_constituent.wt_AA, cds_mutant_constituent.position, cds_mutant_constituent.mut_AA))
            database_ID = plasmid.get_id()
            mutant_ID = cds_mutant.Mutant_ID

        Mutant_plasmid_sequence, CDS_mutant_constituents, mutant_genbank_file = self.generate_mutant_sequence( Plasmid_Feature_ID, mutation_tuples, write_to_file=True)

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






