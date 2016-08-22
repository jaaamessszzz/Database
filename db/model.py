from sqlalchemy import Table, ForeignKey, Column
from sqlalchemy.dialects.mysql import DOUBLE, TINYINT
from sqlalchemy.types import DateTime, Enum, Float, TIMESTAMP, Integer, Text, Unicode, String

# declarative_base instance used for TurboGears integration
from sqlalchemy.ext.declarative import declarative_base
DeclarativeBasePlasmid = declarative_base()


class DBConstants(DeclarativeBasePlasmid):
    __tablename__ = 'DBConstants'

    Parameter = Column(Unicode(256, collation = "utf8_bin"), nullable = False, primary_key = True)
    ParameterType = Column(Unicode(256, collation = "utf8_bin"), nullable = True)
    Value = Column(Text(collation="utf8_bin"), nullable=False)
    ValueType = Column(Enum('Int','String','Float', 'BLOB', 'Timestamp'), nullable = False)


class Primers(DeclarativeBasePlasmid):
    __tablename__ = 'Primers'

    creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Users.ID'), nullable=False, primary_key = True)
    creator_entry_number = Column(Integer, nullable=False, primary_key = True)
    secondary_id = Column(Unicode(100, collation = "utf8_bin"), nullable = False)
    sequence = Column(Unicode(256, collation="utf8_bin"), nullable = False)
    direction = Column(Enum('F','R'), nullable=False)
    description = Column(Text(collation="utf8_bin"), nullable=False)
    TM = Column(Float, nullable=False)
    date = Column(TIMESTAMP, nullable=False)
    GC = Column(Float, nullable=True)
    length = Column(Float, nullable=True)
    template_id = Column(Unicode(200, collation="utf8_bin"), nullable = True)
    template_description = Column(Unicode(200, collation="utf8_bin"), nullable = True)

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
    from sqlalchemy.sql.sqltypes import LargeBinary
    __tablename__ = 'Plasmid'
    creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Users.ID'), nullable=False, primary_key=True)
    creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    plasmid_name = Column(Unicode(100, collation="utf8_bin"), nullable=True, unique=True)
    plasmid_type = Column(Enum("part","cassette","multicassette","other"), nullable=False)
    location = Column(Unicode(250, collation="utf8_bin"), nullable=False)
    description = Column(Text(), nullable=False)
    sequence = Column(Text(), nullable=False)
    ape_file = Column(LargeBinary(), nullable=True)
    status = Column(Enum("designed","verified","abandoned"), nullable=False)
    date = Column(TIMESTAMP, nullable=False)

