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

    def primer_match():
        # Create up the database session
        dbi = DatabaseInterface()
        tsession = dbi.get_session()

        # Create a map from usernames to the database IDs (typically initials)
        user_map = {}
        for u in tsession.query(Users):
            user_map[u.lab_username] = u.ID

        target_sequence = 'ATCAACAACGAGCAGGTCTC'
        plasmid_sequences = tsession.query(Plasmid.plasmid_name, Plasmid.sequence)

        pprint.pprint(plasmid_sequences)

        for name, sequence in plasmid_sequences:
            print name

            if target_sequence in sequence:
                print 'Target sequence was found in %s' %name

Plasmid_Utilities.primer_match()