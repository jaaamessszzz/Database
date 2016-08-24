import os
import sys
import pprint

sys.path.insert(0, '..')

from db.interface import DatabaseInterface
from db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type


class Plasmid_Utilities(object):
    '''
    This is a class with modules for completing tasks involving manipulation of Plasmid sequences.

    ###########################
    # Available Funtionality: #
    ###########################

    * Find available cloning primers for a plasmid of interest
    * Golden Gate Assembly within a modular cloning (MoClo) system - see http://pubs.acs.org/doi/abs/10.1021/sb500366v
    * Verify that Plasmids adhere to MoClo sequence standards

    ###########################
    # Planned Funtionality:   #
    ###########################

    * Generate annotated .ape files based on desired sequence features

    '''
    def __init__(self):
        asdf = 'asdf'

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

    def golden_gate_assembly(self, plasmid_list, restriction_enzyme):
        '''
        This function will take user input of a list of plasmids and attempt to perform a golden gate assembly.

        :param part_plasmid_list: A list of part plasmids that the user wants to assemble into a cassette plasmid
        :param restriction_enzyme: BsaI (cassette assembly) or BsmBI (Multicassette assembly)
        :return: complete plasmid assembly as string
        '''
        # print plasmid_list

        # Create up the database session
        dbi = DatabaseInterface()
        tsession = dbi.get_session()

        # Get sequences and part numbers for listed part plasmids
        part_plasmids_query = tsession.query(Plasmid, Part_Plasmid, Part_Plasmid_Part, Part_Type)\
            .filter(Part_Plasmid.creator == Plasmid.creator)\
            .filter(Part_Plasmid.creator_entry_number == Plasmid.creator_entry_number)\
            .filter(Part_Plasmid.creator == Part_Plasmid_Part.creator)\
            .filter(Part_Plasmid.creator_entry_number == Part_Plasmid_Part.creator_entry_number)\
            .filter(Plasmid.plasmid_name.in_(plasmid_list)) \
            .filter(Part_Plasmid_Part.part_number == Part_Type.part_number)

        if restriction_enzyme == 'BsaI':
            site_F = 'GGTCTC'
            site_R = 'GAGACC'
        if restriction_enzyme == 'BsmBI':
            site_F = 'CGTCTC'
            site_R = 'GAGACG'

        part_list = []
        for plasmid, part_plasmid, part_plasmid_part, part_type in part_plasmids_query:

            print plasmid.plasmid_name, plasmid.creator, plasmid.creator_entry_number, part_plasmid_part.part_number
            print part_type.overhang_5
            print part_type.overhang_3

            assert plasmid.sequence.count(site_F) == 1,\
                'Error: There is more than one forward %s site in %s!' % (restriction_enzyme, plasmid.plasmid_name)
            assert plasmid.sequence.count(site_R) == 1, \
                'Error: There is more than one reverse %s site in %s!' % (restriction_enzyme, plasmid.plasmid_name)

            sequence_upper = plasmid.sequence.upper()
            site_F_position = sequence_upper.find(part_type.overhang_5, sequence_upper.find('GGTCTC'))
            site_R_position = sequence_upper.find(part_type.overhang_3, sequence_upper.find('GAGACC') - 10) + 4

            # Circular permutation so that forward cut site is alway upstream of reverse cut site in linear sequence
            if site_F_position > site_R_position:
                sequence_upper = sequence_upper[site_R_position + 7:] + sequence_upper[:site_R_position + 7]

            part_list.append( sequence_upper[
                              sequence_upper.find(part_type.overhang_5, sequence_upper.find('GGTCTC')):
                              sequence_upper.find(part_type.overhang_3, sequence_upper.find('GAGACC') - 10) + 4
                              ])

        # Start with Part one (cassette assembly) or backbone (multicassette assembly) so that there is consistency in how the plasmids are assembled
        if restriction_enzyme == 'BsaI':
            for part in part_list:
                if part[:4].upper() == 'CCCT':
                    intermediate = part.upper()
                    part_list.remove(part)
                    break
        if restriction_enzyme == 'BsmBI':
            for part in part_list:
                if part[:4].upper() == 'CTGA':
                    intermediate = part.upper()
                    part_list.remove(part)
                    break

        assembly_step = 0
        assembly_steps = len(part_list)

        while assembly_step < assembly_steps:
            for part in part_list:
                # If part 5' and intermediate 3' match...
                if part.upper()[:4] == intermediate.upper()[-4:]:
                    intermediate = intermediate.upper()[:-4] + part.upper()
            assembly_step += 1

        assert intermediate[-4:] == intermediate[:4], 'Error: Incomplete assembly! This assembly does not produce a circular plasmid! :('
        complete_assembly = intermediate[:-4].upper()

        return complete_assembly

    def plasmid_checks(self, plasmid):
        '''
        Verifies that a plasmid conforms to modular cloning sequence standards

        :return:
        '''

