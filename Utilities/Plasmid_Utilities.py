import os
import sys
import pprint

sys.path.insert(0, '..')

from db import model
from db.interface import DatabaseInterface
from db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type


class Plasmid_Utilities(object):
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

    def cassette_assembly(self, part_plasmid_list):
        '''
        This function will take user input of a list of part plasmids and attempt to assemble them into a cassette plasmid.
        Cassette plasmid records should be automatically derived and pushed to the database upon a successful assembly.

        :param part_plasmid_list: A list of part plasmids that the user wants to assemble into a cassette plasmid
        :return:
        '''
        # print part_plasmid_list

        # Create up the database session
        dbi = DatabaseInterface()
        tsession = dbi.get_session()

        # Get sequences and part numbers for listed part plasmids
        part_plasmids_query = tsession.query(Plasmid, Part_Plasmid, Part_Plasmid_Part, Part_Type)\
            .filter(Part_Plasmid.creator == Plasmid.creator)\
            .filter(Part_Plasmid.creator_entry_number == Plasmid.creator_entry_number)\
            .filter(Part_Plasmid.creator == Part_Plasmid_Part.creator)\
            .filter(Part_Plasmid.creator_entry_number == Part_Plasmid_Part.creator_entry_number)\
            .filter(Plasmid.plasmid_name.in_(part_plasmid_list)) \
            .filter(Part_Plasmid_Part.part_number == Part_Type.part_number)



        part_list = []
        for plasmid, part_plasmid, part_plasmid_part, part_type in part_plasmids_query:

            print plasmid.plasmid_name
            print plasmid.creator, plasmid.creator_entry_number
            print part_plasmid_part.part_number
            print part_type.overhang_3
            print part_type.overhang_5

            assert plasmid.sequence.count('GGTCTC') == 1
            assert plasmid.sequence.count('GAGACC') == 1

            sequence_upper = plasmid.sequence.upper()
            part_list.append( plasmid.sequence[
                              sequence_upper.find(part_type.overhang_5, sequence_upper.find('GGTCTC')):
                              sequence_upper.find(part_type.overhang_3, sequence_upper.find('GAGACC') - 10) + 4
                              ])

        for part in part_list:
            if part[:4].upper() == 'CCCT':
                intermediate = part.upper()
                part_list.remove(part)
                break

        while len(part_list) != 0:
            for part in part_list:
                # If part 5' and intermediate 3' match...
                if part.upper()[:4] == intermediate.upper()[-4:]:
                    intermediate = intermediate.upper()[:-4] + part.upper()
                    part_list.remove(part)
                # # If part 3' and intermediate 5' match...
                # elif part.upper()[-4:] == intermediate.upper()[:4]:
                #     intermediate = part.upper() + intermediate.upper()[4:]
                #     part_list.remove(part)
        assert len(part_list) == 0
        complete_assembly = intermediate[:-4].upper()
        pprint.pprint(complete_assembly)

