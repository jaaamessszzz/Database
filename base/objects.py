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