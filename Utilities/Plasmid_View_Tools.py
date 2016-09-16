import sys
import StringIO
import re
import os
from sqlalchemy import and_, or_
from Utilities.Plasmid_Utilities import Plasmid_Utilities, Plasmid_Exception

sys.path.insert(0, '..')

try:
    from db.interface import DatabaseInterface
    from db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature, Plasmid_Feature, Plasmid_File, CDS_Mutant, CDS_Mutant_Constituent
except:
    # nasty hack since we are not packaging things up properly yet for external use (e.g. the website)
    from kprimers.db.interface import DatabaseInterface
    from kprimers.db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature, Plasmid_Feature, Plasmid_File, CDS_Mutant, CDS_Mutant_Constituent

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
        self.plasmid_util = Plasmid_Utilities(tsession or self.dbi.get_session())

        # Get User ID
        self.username = username
        if not self.username:
            self.username = os.getlogin() # The webserver cannot call this function so it must instead pass username in as an argument
        self.user_ID = self.tsession.query(Users).filter(Users.lab_username == self.username).one().ID

    def get_plasmid_indicies(self, current_plasmid):
        """
        Generates two lists that contain starting and ending indicies for parts and features within a given plasmid. \
        These are used by the AngularPlasmid javascript library to display an annotated plasmid map.

        current_plasmid: tuple of (creator, creator_entry_number) for current plasmid

        :return:
        feature_indicies_list: list of name, tuple with start and end for each feature
        part_indicies_list: list of name, tuple with start and end for each part
        """

        current_plasmid_sequence = self.tsession.query(Plasmid).filter(and_(Plasmid.creator == current_plasmid[0],
                                                                            Plasmid.creator_entry_number == current_plasmid[1])
                                                                       ).one()

        plasmid_features = self.tsession.query(Plasmid_Feature, Feature, Feature_Type)\
            .filter(and_(Plasmid_Feature.creator == current_plasmid_sequence.creator,
                         Plasmid_Feature.creator_entry_number == current_plasmid_sequence.creator_entry_number,
                         Plasmid_Feature.feature_name == Feature.Feature_name,
                         Feature_Type.Feature_type == Feature.Feature_type
                         )
                    )

        feature_indicies_list = [] # List of [feature_name, (starting_index, ending_index)]

        # todo: add all relevant information to feature index list (only name and indicies so far)
        for plasmid_feature, feature, feature_type in plasmid_features:
            for instance in re.finditer(feature.Feature_sequence.upper().strip(), current_plasmid_sequence.sequence):
                if feature.Feature_type == 'Restrxn Type II':
                    feature_indicies_list.append([plasmid_feature.feature_name + ' F', (instance.start(), instance.end())])
                    feature_indicies_list.append([plasmid_feature.feature_name + ' Overhang', (instance.start() + 7, instance.start() + 11)])
                else:
                    feature_indicies_list.append([plasmid_feature.feature_name, (instance.start(), instance.end())])

            for instance in re.finditer(self.plasmid_util.reverse_complement(feature.Feature_sequence.upper().strip()), current_plasmid_sequence.sequence):
                if feature.Feature_type == 'Restrxn Type II':
                    feature_indicies_list.append([plasmid_feature.feature_name + ' R', (instance.start(), instance.end())])
                    feature_indicies_list.append([plasmid_feature.feature_name + ' Overhang', (instance.start() - 5, instance.start() - 1)])
                else:
                    feature_indicies_list.append([plasmid_feature.feature_name, (instance.start(), instance.end())])

        # todo: add all relevant information to part index list (only name and indicies so far)
        if current_plasmid_sequence.plasmid_type.lower() == 'cassette':
            # List of [part_name, (starting_index, ending_index)]
            part_indicies_list = []

            cassette_parts = self.tsession.query(Cassette_Assembly, Plasmid)\
                .filter(and_(current_plasmid_sequence.creator == Cassette_Assembly.Cassette_creator,
                             current_plasmid_sequence.creator_entry_number == Cassette_Assembly.Cassette_creator_entry_number,
                             Plasmid.creator == Cassette_Assembly.Part_creator,
                             Plasmid.creator_entry_number == Cassette_Assembly.Part_creator_entry_number))

            for cassette_assembly, plasmid in cassette_parts:
                # Do the restriction digest with BsaI and get indices for start/end of part minus overhangs
                part_sequence = self.plasmid_util.fetch_cassette_parts([(plasmid.creator, plasmid.creator_entry_number)])[4:-4]

                for part in part_sequence:
                    for instance in re.finditer(part.upper().strip(), current_plasmid_sequence.sequence):
                        part_indicies_list.append([plasmid.description, (instance.start(), instance.end())])

        return feature_indicies_list, part_indicies_list