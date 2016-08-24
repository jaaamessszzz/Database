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

    creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Users.ID'), nullable=False, primary_key=True)
    creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    plasmid_name = Column(Unicode(100, collation="utf8_bin"), nullable=True, unique=True)
    plasmid_type = Column(Enum("part","cassette","multicassette","other"), nullable=False)
    location = Column(Unicode(250, collation="utf8_bin"), nullable=False)
    description = Column(Text(), nullable=False)
    sequence = Column(Text(), nullable=False)
    status = Column(Enum("designed","verified","abandoned"), nullable=False)
    date = Column(TIMESTAMP, nullable=False)

    # def __repr__(self):
    #     return 'Internal numbering: {0} {1}\n' \
    #           'plasmid_name: {2}'.format(self.creator, self.creator_entry_number, self.plasmid_name)


######

class Cassette_Assembly(DeclarativeBasePlasmid):
    __tablename__ = 'Cassette_Assembly'

    Cassette_creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    Cassette_creator = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    Part_number = Column(String(2), nullable=False, primary_key=True)
    Part_creator = Column(Unicode(5, collation="utf8_bin"), nullable=False)
    Part_creator_entry_number = Column(Integer, nullable=False)


class Cassette_Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Cassette_Plasmid'

    creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    creator = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    left_connector = Column(Integer, nullable=True)
    right_connector = Column(Integer, nullable=True)


class Feature(DeclarativeBasePlasmid):
    __tablename__ = 'Feature'

    Feature_name = Column(Unicode(200, collation="utf8_bin"), nullable=False, primary_key=True)
    MD5_hash = Column(String(64), nullable=True)
    Feature_type = Column(Unicode(200, collation="utf8_bin"), nullable=False)
    Feature_sequence = Column(Text, nullable=True)


class Feature_Type(DeclarativeBasePlasmid):
    __tablename__ = 'Feature_Type'

    Feature_type = Column(Unicode(200, collation="utf8_bin"), nullable=False, primary_key=True)


class Left_Cassette_Connector(DeclarativeBasePlasmid):
    __tablename__ = 'Left_Cassette_Connector'

    connector = Column(Integer, nullable=False, primary_key=True)
    overhang = Column(String(4), nullable=True)


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

    creator_entry_number = Column(Integer, ForeignKey('Plasmid.creator_entry_number'), nullable=False, primary_key=True)
    creator = Column(Unicode(5), ForeignKey('Plasmid.creator'), nullable=False, primary_key=True)
    resistance = Column(Unicode(100, collation="utf8_bin"), nullable=False)

    # ID = relationship("Plasmid", primaryjoin='and_(Plasmid.creator==Part_Plasmid.creator, Plasmid.creator_entry_number==Part_Plasmid.creator_entry_number)')

    # def __repr__(self):
    #     return 'Internal numbering: {0} {1}'.format(self.creator, self.creator_entry_number)

class Part_Plasmid_Part(DeclarativeBasePlasmid):
    __tablename__ = 'Part_Plasmid_Part'

    creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    creator = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    part_number = Column(String(4), nullable=False, primary_key=True)

    # def __repr__(self):
    #     return 'Internal numbering: {0} {1}\nPart number: {2}'.format(self.creator, self.creator_entry_number, self.part_number)

class Part_Type(DeclarativeBasePlasmid):
    __tablename__ = 'Part_Type'

    part_number = Column(Unicode(4), nullable=False, primary_key=True)
    overhang_5 = Column(Unicode(4), nullable=False)
    overhang_3 = Column(Unicode(4), nullable=False)


class Plasmid_Feature(DeclarativeBasePlasmid):
    __tablename__ = 'Plasmid_Feature'

    ID = Column(Integer, nullable=False, primary_key=True)
    creator_entry_number = Column(Integer, nullable=False)
    creator = Column(Unicode(5, collation="utf8_bin"), nullable=False)
    feature_name = Column(Unicode(200, collation="utf8_bin"), nullable=False)


class Right_Cassette_Connector(DeclarativeBasePlasmid):
    __tablename__ = 'Right_Cassette_Connector'

    connector = Column(Integer, nullable=False, primary_key=True)
    overhang = Column(String(4), nullable=True)


