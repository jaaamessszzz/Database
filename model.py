from sqlalchemy import Table, ForeignKey, Column
from sqlalchemy.dialects.mysql import DOUBLE, TINYINT
from sqlalchemy.types import DateTime, Enum, Float, TIMESTAMP, Integer, Text, Unicode

# declarative_base instance used for TurboGears integration
from sqlalchemy.ext.declarative import declarative_base
DeclarativeBasePlasmid = declarative_base()


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



class Resistance(DeclarativeBasePlasmid):
    __tablename__ = 'Resistance'

    resistance = Column(Unicode(100, collation="utf8_bin"), nullable=False, primary_key=True)


class Users(DeclarativeBasePlasmid):
    __tablename__ = 'Users'

    ID = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    github_username = Column(Unicode(26, collation="utf8_bin"), nullable=False)
    lab_username = Column(Unicode(50, collation="utf8_bin"), nullable = False)
