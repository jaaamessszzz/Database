import sys
import StringIO
import re
import os
import numpy as np
from sqlalchemy import and_, or_

sys.path.insert(0, '..')

try:
    from db.interface import DatabaseInterface
    from db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature, Plasmid_Feature, Plasmid_File, CDS_Mutant, CDS_Mutant_Constituent
    from Utilities.Plasmid_Utilities import Plasmid_Utilities, Plasmid_Exception
except:
    # nasty hack since we are not packaging things up properly yet for external use (e.g. the website)
    from kprimers.db.interface import DatabaseInterface
    from kprimers.db.model import Users, Plasmid, Primers, Part_Plasmid, Part_Plasmid_Part, Part_Type, Cassette_Assembly, Cassette_Plasmid, Cassette_Connector, Feature_Type, Feature, Plasmid_Feature, Plasmid_File, CDS_Mutant, CDS_Mutant_Constituent
    from kprimers.Utilities.Plasmid_Utilities import Plasmid_Utilities, Plasmid_Exception


class PrimerGenerationException(Exception): pass


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
        self.plasmid_util = Plasmid_Utilities(tsession or self.dbi.get_session(), username = self.username)


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
            for instance in re.finditer(feature.Feature_sequence.upper().strip(), current_plasmid_sequence.sequence.upper()):
                if feature.Feature_type == 'Restrxn Type II':
                    feature_indicies_list.append([plasmid_feature.feature_name + ' F', (instance.start(), instance.end())])
                    feature_indicies_list.append([plasmid_feature.feature_name + ' Overhang', (instance.start() + 7, instance.start() + 11)])
                else:
                    feature_indicies_list.append([plasmid_feature.feature_name, (instance.start(), instance.end())])

            for instance in re.finditer(self.plasmid_util.reverse_complement(feature.Feature_sequence.upper().strip()), current_plasmid_sequence.sequence.upper()):
                if feature.Feature_type == 'Restrxn Type II':
                    feature_indicies_list.append([plasmid_feature.feature_name + ' R', (instance.start(), instance.end()), feature_type.color])
                    feature_indicies_list.append([plasmid_feature.feature_name + ' Overhang', (instance.start() - 5, instance.start() - 1), feature_type.color])
                else:
                    feature_indicies_list.append([plasmid_feature.feature_name, (instance.start(), instance.end()), feature_type.color])

        # todo: add all relevant information to part index list (only name and indicies so far)
        part_indicies_list = []
        if current_plasmid_sequence.plasmid_type.lower() == 'cassette':
            # List of [part_name, (starting_index, ending_index)]

            cassette_parts = self.tsession.query(Cassette_Assembly, Plasmid)\
                .filter(and_(current_plasmid_sequence.creator == Cassette_Assembly.Cassette_creator,
                             current_plasmid_sequence.creator_entry_number == Cassette_Assembly.Cassette_creator_entry_number,
                             Plasmid.creator == Cassette_Assembly.Part_creator,
                             Plasmid.creator_entry_number == Cassette_Assembly.Part_creator_entry_number))

            for cassette_assembly, plasmid in cassette_parts:
                # Do the restriction digest with BsaI and get indices for start/end of part minus overhangs
                part_sequence = self.plasmid_util.fetch_cassette_parts([(plasmid.creator, plasmid.creator_entry_number)])[4:-4]

                for part in part_sequence:
                    for instance in re.finditer(part.upper().strip(), current_plasmid_sequence.sequence.upper()):
                        part_indicies_list.append([plasmid.description, (instance.start(), instance.end())])

        return feature_indicies_list, part_indicies_list


    def generate_primers(self, target_sequence, Target_TM, primer = 50, Na = 50, K = None, Mg = None, dNTPs = None, Tris = None, left_arm = None, right_arm = None):
        '''
        Function for generating primers with SequenceViewer + User inputs and depositing designed primers into the
        Primer table.

        Tried doing this by hand but discovered that BioPython already has a function for this...

        If so inclined, here are the references I used:
        d_H and d_S values from Improved thermodynamic parameters and helix initiation factor to predict stability of
        DNA duplexes. Naoki Sugimoto, Shu-ich Nakano, Mari Yoneyama and Kei-ich Honda. Nucleic Acids Research, 1996,
        Vol. 24, No. 22 4501-4505
        Equations and stuff from A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor
        thermodynamics. JOHN SANTALUCIA, JR. Proceedings of the National Academy of Sciences. Vol. 95, pp 1460-1465,
        February 1998

        :param target_sequence: Input target sequence from SequenceViewer
        :param TM: Target melting temperature in Celcius
        :param salt: salt concentration in Molar
        :return: Forward and Reverse Primers with associated TMs
        '''

        from Bio.SeqUtils import MeltingTemp as mt
        from Bio.Seq import Seq

        primer_length_F = 1
        primer_length_R = 1

        TM_F = 0
        TM_R = 0

        target_upper = target_sequence.upper()
        for char in target_upper:
            if char.upper() not in 'ATCG':
                raise Plasmid_Exception('Non-ATCG character in the target sequence!')

        # Generate Forward Primer
        target_primer_F, TM_F_previous = None, None
        while TM_F < Target_TM:
            TM_F = mt.Tm_NN(Seq(target_upper[:primer_length_F]),
                            nn_table = mt.DNA_NN2,
                            dnac1 = primer / 2,  # nM Primers / 2
                            dnac2 = primer / 2,  # nM Primers / 2
                            selfcomp = False,
                            Na = Na,     # mM
                            K = K or 0,        # mM
                            Tris = Tris or 0,     # mM
                            Mg = Mg or 0,       # mM
                            dNTPs = dNTPs or 0,
                            saltcorr = 5
                            )
            if TM_F < Target_TM:
                primer_length_F += 1
                target_primer_F = target_upper[:primer_length_F]
            if TM_F_previous != None and TM_F_previous >= TM_F:
                # Strong normalisation check - otherwise this will loop infinitely (possible with bad input)
                target_primer_F = None
                break
            TM_F_previous = TM_F

        if target_primer_F == None:
            raise PrimerGenerationException()

        # Generate Reverse Primer
        target_primer_R, TM_R_previous = None, None
        while TM_R < Target_TM:
            TM_R = mt.Tm_NN(Seq(self.plasmid_util.reverse_complement(target_upper[-primer_length_R:])),
                            nn_table=mt.DNA_NN2,
                            dnac1=primer / 2,  # nM Primers / 2
                            dnac2=primer / 2,  # nM Primers / 2
                            selfcomp=False,
                            Na=Na,  # mM
                            K=K or 0,  # mM
                            Tris=Tris or 0,  # mM
                            Mg=Mg or 0,  # mM
                            dNTPs=dNTPs or 0,
                            saltcorr=5
                            )
            if TM_R < Target_TM:
                primer_length_R += 1
                target_primer_R = self.plasmid_util.reverse_complement(target_upper[-primer_length_R:])
            if TM_R_previous != None and TM_R_previous >= TM_R:
                # Strong normalisation check - otherwise this will loop infinitely (possible with bad input)
                target_primer_R = None
                break
            TM_R_previous = TM_R

        if target_primer_R == None:
            raise PrimerGenerationException()

        if left_arm and right_arm:
            for char in left_arm:
                if char.upper() not in 'ATCG':
                    raise Plasmid_Exception('Non-ATCG character in the primer extensions!')
            target_primer_F = left_arm.upper() + target_primer_F
            target_primer_R = self.plasmid_util.reverse_complement(right_arm).upper() + target_primer_R

            return target_primer_F, target_primer_R
        elif left_arm == None and right_arm == None:
            return target_primer_F, target_primer_R
        else:
            raise Plasmid_Exception('You need to add both the left and right extensions for this to work!')



