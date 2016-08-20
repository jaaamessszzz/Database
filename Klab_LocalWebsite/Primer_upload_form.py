import os
import sys
import pprint
import pandas as pd

import sys
sys.path.insert(0, '..')

from model import Primers

test_csv_path = 'biosensor_primers_v2.csv'
delimiter = '\t'

print os.getcwd()

# Import data
user_input = pd.read_csv(test_csv_path, delimiter='\t')

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


# Parse and add to model
print(user_input)
for index, row in user_input.iterrows():

    d = {}
    for header in ['creator', 'secondary_id', 'sequence', 'direction', 'description', 'TM', 'GC', 'length', 'template_id', 'template_description']:
        d[header] = row[header]
    db_record_object = Primers(**d)

    pprint.pprint(vars(Primers))




