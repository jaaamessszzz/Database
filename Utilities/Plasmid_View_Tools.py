import sys
import StringIO
import re
import os
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

class Plasmid_View_Tools(object):
    '''
    Tools for viewing and manipulating plasmids in Plasmid View!
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

    def get_plasmid_indicies(self, current_plasmid):
        """
        Generates an object (to be determined!) that contains starting and ending indicies for parts and features \
        contained within a given plasmid. These are used by the AngularPlasmid javascript library to display an \
        annotated plasmid map.

        current_plasmid: tuple of (creator, creator_entry_number) for current plasmid

        :return:
        """

        current_plasmid_sequence = self.tsession.query(Plasmid).filter(and_(Plasmid.creator == current_plasmid[0],
                                                                            Plasmid.creator_entry_number == current_plasmid[1]
                                                                            )
                                                                       ).one()

        print current_plasmid_sequence