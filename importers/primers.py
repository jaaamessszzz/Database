import os
import sys
import pprint
import traceback
import datetime

import pandas as pd

from klab import colortext

import sys
sys.path.insert(0, '..')

from db import model
from db.interface import DatabaseInterface
from db.model import DBConstants, Primers, Users


def parse(primers_file, user_id, errors = []):
    print(primers_file)
    #tsession = dbi.get_session(new_session = True)

    # Import data
    df = pd.read_csv(primers_file, delimiter='\t')
    colortext.warning(str(df))


    required_fields = ['sequence', 'direction', 'description', 'template_id']
    optional_fields = ['secondary_id', 'TM', 'GC', 'length', 'template_description']
    derived_or_generated_fields = ['creator', 'creator_entry_number', 'date']
    float_fields = ['TM', 'GC', 'length']

    # Parse and add to model
    errors = []
    print(df.columns.values)
    for f in required_fields:
        if f not in df.columns.values:
            errors.append('Missing required field {0} in {1}.'.format(f, primers_file))
    for f in derived_or_generated_fields:
        if f in df.columns.values:
            errors.append('Found derived field {0} in {1}. Since these values are derived during import into the database, storing them explicitly in the file could lead to inconsistencies.'.format(f, primers_file))

    if errors:
        raise Exception('Failed to parse the TSV file.')

    # Set up derived fields. creator_entry_number is created automatically by a database trigger upon import

    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    for index, row in df.iterrows():

        d = dict(
            creator = user_id,
            date = timestamp,
        )

        for header in required_fields:
            d[header] = row[header]
        for header in optional_fields:
            if header in row:
                d[header] = row[header]
            else:
                d[header] = None
        for header in float_fields:
            if d[header] != None:
                d[header] = float(d[header])

        colortext.pcyan('Adding this record:')
        pprint.pprint(d)

        db_record_object = Primers(**d)

        colortext.warning('String representation:')
        print(db_record_object)
        print('')


def main():

    # Create up the database session
    dbi = DatabaseInterface()
    tsession = dbi.get_session()

    # Create a map from usernames to the database IDs (typically initials)
    user_map = {}
    for u in tsession.query(Users):
        user_map[u.lab_username] = u.ID

    # Read the import path from the database
    colortext.message('\nPrimers import script')
    colortext.pcyan('Database admin contacts: {0}'.format(', '.join(dbi.get_admin_contacts())))
    colortext.warning('Registered users: {0}\n'.format(', '.join(   ['{0} ({1})'.format(v, k) for k, v in sorted(user_map.iteritems(), key = lambda x: x[1])])))
    
    errors = []
    import_path = tsession.query(DBConstants).filter(DBConstants.Parameter == u'import_path').one().Value
    import_path_folders = sorted([d for d in os.listdir(import_path) if os.path.isdir(os.path.join(import_path,d))])
    for ipf in import_path_folders:
        if ipf in user_map:
            user_folder = os.path.join(import_path, ipf)
            user_id = user_map[ipf]
            primers_file = os.path.join(user_folder, 'primers.tsv')
            if os.path.exists(primers_file):
                try:
                    parse(primers_file, user_id, errors)
                except Exception, e:
                    colortext.warning('Error: {0}\n{1}'.format(str(e), traceback.format_exc()))

    if errors:
        pass
        # todo: email the plasdmins


if __name__ == '__main__':
    main()





