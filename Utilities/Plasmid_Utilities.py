import os
import sys
import pprint

sys.path.insert(0, '..')

from db import model
from db.interface import DatabaseInterface
from db.model import Users, Plasmid


class Plasmid_Utilities(object):
    def __init__(self):
        asdf = 'asdf'

    def reverse_complement(self, input_sequence):
        base_pair = {'A': 'T',
                     'T': 'A',
                     'C': 'G',
                     'G': 'C'
                     }
        return ''.join(reversed([base_pair[base] for base in input_sequence]))

    def primer_match(self, target_sequence):
        # Create up the database session
        dbi = DatabaseInterface()
        tsession = dbi.get_session()

        # Create a map from usernames to the database IDs (typically initials)
        user_map = {}
        for u in tsession.query(Users):
            user_map[u.lab_username] = u.ID

        plasmid_sequences = tsession.query(Plasmid.plasmid_name, Plasmid.sequence)

        target_reverse_complement = self.reverse_complement(target_sequence)

        for name, sequence in plasmid_sequences:
            print name
            if target_sequence in sequence:
                print 'Target sequence (F) was found in %s' % name
            if target_reverse_complement in sequence:
                print 'Target sequence (R) was found in %s' % name
