"""
This file was written by James Lucas and Shane O'Connor as part of the Kortemme Lab Plasmid Database.
"""
__author__ = "James Lucas and Shane O'Connor"
__copyright__ = ""
__credits__ = ["James Lucas", "Shane O'Connor"]
__license__ = "CC Attribution-NonCommercial 4.0 International"
__version__ = "0.1"
__maintainer__ = "James Lucas"
__email__ = "james.lucas@berkeley.edu"
__status__ = "Development"

class user_input(object):
    def __init__(self, creator, input_sequences, assembly_type, UID, location, part_number = None,
                 part_description = None):
        self.creator = creator
        self.input_sequences = input_sequences
        self.assembly_type = assembly_type
        self.UID = UID
        self.location = location
        self.part_number = part_number
        self.part_description = part_description