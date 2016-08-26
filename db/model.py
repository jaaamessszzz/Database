import pprint

import numpy

from sqlalchemy import Table, ForeignKey, Column
from sqlalchemy.dialects.mysql import DOUBLE, TINYINT, LONGBLOB
from sqlalchemy.types import DateTime, Enum, Float, TIMESTAMP, Integer, Text, Unicode, String
from sqlalchemy.orm import relationship

# declarative_base instance used for TurboGears integration
from sqlalchemy.ext.declarative import declarative_base
DeclarativeBasePlasmid = declarative_base()

from klab import colortext


class DBConstants(DeclarativeBasePlasmid):
    __tablename__ = 'DBConstants'

    Parameter = Column(Unicode(256, collation = "utf8_bin"), nullable = False, primary_key = True)
    ParameterType = Column(Unicode(256, collation = "utf8_bin"), nullable = True)
    Value = Column(Text(collation="utf8_bin"), nullable=False)
    ValueType = Column(Enum('Int','String','Float', 'BLOB', 'Timestamp'), nullable = False)


class Primers(DeclarativeBasePlasmid):
    __tablename__ = 'Primers'

    _required_fields = ['sequence', 'direction', 'description', 'TM']
    _optional_fields = ['secondary_id', 'GC', 'length', 'template_id', 'template_description']
    _derived_or_generated_fields = ['creator', 'creator_entry_number', 'date']
    _float_fields = ['GC', 'TM']
    _int_fields = ['length']

    creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Users.ID'), nullable=False, primary_key = True)
    creator_entry_number = Column(Integer, nullable=False, primary_key = True)
    secondary_id = Column(Unicode(100, collation = "utf8_bin"), nullable = False)
    sequence = Column(Unicode(256, collation="utf8_bin"), nullable = False)
    direction = Column(Enum('F','R'), nullable=False)
    description = Column(Text(collation="utf8_bin"), nullable=False)
    TM = Column(Integer, nullable=False)
    date = Column(TIMESTAMP, nullable=False)
    GC = Column(Float, nullable=True)
    length = Column(Float, nullable=True)
    template_id = Column(Unicode(200, collation="utf8_bin"), nullable = True)
    template_description = Column(Unicode(200, collation="utf8_bin"), nullable = True)


    @staticmethod
    def add(tsession, d, silent = True):

        try:
            for f in Primers._optional_fields:
                if f not in d:
                    d[f] = None
            for f in Primers._float_fields:
                if d[f] != None:
                    d[f] = float(d[f])
            for f in Primers._int_fields:
                if d[f] != None:
                    d[f] = int(d[f])
            for k, v in d.iteritems():
                if type(v) == float and numpy.isnan(v):
                    d[k] = None
            db_record_object = Primers(**d)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()
        except:
            raise


    def __repr__(self):
        description = self.description or ''
        if description:
            description = '\n  Description: {0}'.format(description)
        GC = self.GC or ''
        if GC:
            GC = ', {0}'.format(GC)
        plength = self.length or ''
        if plength:
            plength = ', {0}'.format(plength)
        if self.creator_entry_number:
            # This is added on database insertion so will not be filled in for transient objects
            s = 'Primer o{0}{1:04d} ({2}): {3}, {6}C, {7}{8} \n  Sequence: {4}{5}'.format(self.creator, self.creator_entry_number, self.secondary_id,
                                                                    self.direction, self.sequence, description, self.TM, self.date, GC, plength)
        else:
            s = 'Primer o{0} ({2}): {3}, {6}C, {7}{8} \n  Sequence: {4}{5}'.format(self.creator, None, self.secondary_id,
                                                                    self.direction, self.sequence, description, self.TM, self.date, GC, plength)
        if self.template_id:
            s += '\n  Template #: {0}'.format(self.template_id)
        if self.template_description:
            s += '\n  Template: {0}'.format(self.template_description)
        return s


class Resistance(DeclarativeBasePlasmid):
    __tablename__ = 'Resistance'

    resistance = Column(Unicode(100, collation="utf8_bin"), nullable=False, primary_key=True)


class Users(DeclarativeBasePlasmid):
    __tablename__ = 'Users'

    ID = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    github_username = Column(Unicode(26, collation="utf8_bin"), nullable=False)
    lab_username = Column(Unicode(50, collation="utf8_bin"), nullable = False)
    color = Column(Unicode(6, collation="utf8_bin"), nullable = False)

# TEMP for testing Find_Primers.py
class Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Plasmid'

    _required_fields = ['sequence', 'description', 'location', 'plasmid type', 'status']
    _optional_fields = ['plasmid_name']
    _derived_or_generated_fields = ['creator', 'creator_entry_number', 'date']

    creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Users.ID'), nullable=False, primary_key=True)
    creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    plasmid_name = Column(Unicode(100, collation="utf8_bin"), nullable=True, unique=True)
    plasmid_type = Column(Enum("part","cassette","multicassette","other"), nullable=False)
    location = Column(Unicode(250, collation="utf8_bin"), nullable=False)
    description = Column(Text(), nullable=False)
    sequence = Column(Text(), nullable=False)
    status = Column(Enum("designed","verified","abandoned"), nullable=False)
    date = Column(TIMESTAMP, nullable=False)

    @staticmethod
    def add(tsession, input_dict, silent=True):

        try:
            for f in Plasmid._optional_fields:
                if f not in input_dict:
                    input_dict[f] = None

            pprint.pprint(input_dict)

            db_record_object = Plasmid(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()
            return db_record_object
        except:
            raise

    # def __repr__(self):
    #     return 'Internal numbering: {0} {1}\n' \
    #           'plasmid_name: {2}'.format(self.creator, self.creator_entry_number, self.plasmid_name)

######

class Cassette_Assembly(DeclarativeBasePlasmid):
    __tablename__ = 'Cassette_Assembly'

    _required_fields = ['Cassette_creator_entry_number', 'Cassette_creator', 'Part_number', 'Part_creator', 'Part_creator_entry_number']

    Cassette_creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    Cassette_creator = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    Part_number = Column(String(2), nullable=False, primary_key=True)
    Part_creator = Column(Unicode(5, collation="utf8_bin"), nullable=False)
    Part_creator_entry_number = Column(Integer, nullable=False)

    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Cassette_Assembly(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()
        except:
            raise

class Cassette_Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Cassette_Plasmid'

    _required_fields = ['left_connector_part', 'left_overhang', 'right_connector_part', 'right_overhang type']
    _derived_or_generated_fields = ['creator', 'creator_entry_number']

    creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    creator = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    left_connector_part = Column(Enum('1','5'), nullable=False, primary_key=True)
    left_overhang = Column(Unicode(4, collation="utf8_bin"), nullable=False, primary_key=True)
    right_connector_part = Column(Enum('1','5'), nullable=False, primary_key=True)
    right_overhang = Column(Unicode(4, collation="utf8_bin"), nullable=False, primary_key=True)

    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Cassette_Plasmid(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()
        except:
            raise

class Feature(DeclarativeBasePlasmid):
    __tablename__ = 'Feature'

    Feature_name = Column(Unicode(200, collation="utf8_bin"), nullable=False, primary_key=True)
    MD5_hash = Column(String(64), nullable=True)
    Feature_type = Column(Unicode(200, collation="utf8_bin"), nullable=False)
    Feature_sequence = Column(Text(collation="utf8_bin"), nullable=True)
    description = Column(Text(collation="utf8_bin"))


class Feature_Type(DeclarativeBasePlasmid):
    __tablename__ = 'Feature_Type'

    Feature_type = Column(Unicode(200, collation="utf8_bin"), nullable=False, primary_key=True)
    color = Column(Unicode(7, collation="utf8_bin"), nullable=False)


class Multicassette_Assembly(DeclarativeBasePlasmid):
    __tablename__ = 'Multicassette_Assembly'

    Multicassette_creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    Multicassette_creator = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    Left_connector = Column(Integer, nullable=False, primary_key=True)
    Right_connector = Column(Integer, nullable=False, primary_key=True)
    Cassette_creator_entry_number = Column(Integer, nullable=False)
    Cassette_creator = Column(Unicode(5, collation="utf8_bin"), nullable=False)


class Multicassette_Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Multicassette_Plasmid'

    creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    creator = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    resistance = Column(Unicode(100), nullable=False)
    Origin = Column(Unicode(100, collation="utf8_bin"), nullable=False)


class Origin(DeclarativeBasePlasmid):
    __tablename__ = 'Origin'

    Origin = Column(Unicode(100, collation="utf8_bin"), nullable=False, primary_key=True)


class Part_Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Part_Plasmid'

    _required_fields = ['creator', 'creator_entry_number', 'resistance']

    creator_entry_number = Column(Integer, ForeignKey('Plasmid.creator_entry_number'), nullable=False, primary_key=True)
    creator = Column(Unicode(5), ForeignKey('Plasmid.creator'), nullable=False, primary_key=True)
    resistance = Column(Unicode(100, collation="utf8_bin"), nullable=False)

    # ID = relationship("Plasmid", primaryjoin='and_(Plasmid.creator==Part_Plasmid.creator, Plasmid.creator_entry_number==Part_Plasmid.creator_entry_number)')

    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Part_Plasmid(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()
        except:
            raise

    def __repr__(self):
        return 'Internal numbering: {0} {1}'.format(self.creator, self.creator_entry_number)

class Part_Plasmid_Part(DeclarativeBasePlasmid):
    __tablename__ = 'Part_Plasmid_Part'

    _required_fields = ['creator', 'creator_entry_number', 'part_number']

    creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    creator = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    part_number = Column(String(4), nullable=False, primary_key=True)

    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Part_Plasmid_Part(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()
        except:
            raise

    def __repr__(self):
        return 'Internal numbering: {0} {1}\nPart number: {2}'.format(self.creator, self.creator_entry_number, self.part_number)

class Part_Type(DeclarativeBasePlasmid):
    __tablename__ = 'Part_Type'

    part_number = Column(Unicode(4), nullable=False, primary_key=True)
    overhang_5 = Column(Unicode(4), nullable=False)
    overhang_3 = Column(Unicode(4), nullable=False)


class Plasmid_Feature(DeclarativeBasePlasmid):
    __tablename__ = 'Plasmid_Feature'

    _required_fields = ['creator', 'creator_entry_number', 'feature_name']
    _derived_or_generated_fields = ['ID']

    ID = Column(Integer, nullable=False, primary_key=True)
    creator_entry_number = Column(Integer, nullable=False)
    creator = Column(Unicode(5, collation="utf8_bin"), nullable=False)
    feature_name = Column(Unicode(200, collation="utf8_bin"), nullable=False)

    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Plasmid_Feature(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()
        except:
            raise

class Cassette_Connector(DeclarativeBasePlasmid):
    __tablename__ = 'Cassette_Connector'

    connector_ID = Column(Unicode(2, collation="utf8_bin"), nullable=False, unique=True)
    creator_entry_number = Column(Integer, ForeignKey('Part_Plasmid_Part.creator_entry_number'), nullable=False)
    creator = Column(Unicode(5), ForeignKey('Part_Plasmid_Part.creator'), nullable=False)
    connector_part = Column(Enum('1','5'), nullable=False, primary_key=True)
    overhang = Column(Unicode(4, collation="utf8_bin"), nullable=False, primary_key=True)


class File_Type(DeclarativeBasePlasmid):
    __tablename__ = 'File_Type'

    file_type = Column(Unicode(100), nullable=False, primary_key=True)
    file_extension = Column(String(16), nullable = True)

class Plasmid_File(DeclarativeBasePlasmid):
    __tablename__ = 'Plasmid_File'

    _required_fields = ['creator', 'creator_entry_number', 'file_name', 'file_type', 'Description', 'File']

    ID = Column(Integer, nullable=False, primary_key=True)
    creator_entry_number = Column(Integer, ForeignKey('Part_Plasmid_Part.creator_entry_number'), nullable=False)
    creator = Column(Unicode(5), ForeignKey('Part_Plasmid_Part.creator'), nullable=False)
    file_name = Column(Unicode(100, collation="utf8_bin"), nullable=False)
    file_type = Column(Unicode(100), nullable=False)
    Description = Column( Text(), nullable=False)
    File = Column(LONGBLOB(), nullable=False)

    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Plasmid_File(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()
        except:
            raise
