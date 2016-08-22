import os
import sys
sys.path.insert(0, '..')

from db import model
from db.interface import DatabaseInterface
from db.model import DBConstants, Primers, Users, Plasmid

def main():
    # Create up the database session
    dbi = DatabaseInterface()
    tsession = dbi.get_session()

    # Create a map from usernames to the database IDs (typically initials)
    user_map = {}
    for u in tsession.query(Users):
        user_map[u.lab_username] = u.ID

    target_sequence = 'ATCAACAACGAGCAGGTCTC'
    plasmid_sequences = tsession.query(Plasmid.plasmid_name, Plasmid.sequence)
    for name, sequence in plasmid_sequences:
        print (name, sequence)

if __name__ == '__main__':
    main()