import os
import sys
import pprint
import pandas as pd

test_csv_path = 'biosensor_primers_v2.csv'
delimiter = ','

print os.getcwd()

# Import data
user_input = pd.read_csv(test_csv_path, delimiter=delimiter)

# Model
# class Primers(DeclarativeBasePlasmid):
#     __tablename__ = 'Primers'
#
#     name = Column(Unicode(100), nullable=False, primary_key=True)
#     creator = Column(Unicode(50), nullable=False)
#     direction = Column(Enum('F','R'), nullable=False)
#     description = Column(Text, nullable=False)
#     TM = Column(Float, nullable=False)
#     date = Column(TIMESTAMP, nullable=False)
#     GC = Column(Float, nullable=True)
#     length = Column(Float, nullable=True)
#     template_id = Column(Unicode(200), nullable = True)
#     template_description = Column(Unicode(200), nullable = True)

# Pseudo-model
class Primers():
    def __init__(self):
        self.name = None
        self.creator = None
        self.sequence = None
        self.direction = None
        self.description = None
        self.TM = None
        self.date = None
        self.GC = None
        self.length = None
        self.template_id = None
        self.template_description = None

# Parse and add to model
for index, row in user_input.iterrows():
    Primers.name = row[0]
    Primers.creator = row[1]
    Primers.sequence = row[2]
    Primers.direction = row[3]
    Primers.description = row[4]
    Primers.TM = row[5]
    Primers.date = row[6]
    Primers.GC = row[7]
    Primers.length = row[8]
    Primers.template_id = row[9]
    Primers.template_description = row[10]
    pprint.pprint(vars(Primers))

# from sqlalchemy.orm import sessionmaker
# Session = sessionmaker(bind=engine)
# session = Session()
#
# plasmid = Plasmids(name='pJL001', user='James', sequence='ATCGATCGATCGATCGATCG')
#
# session.add(plasmid)


