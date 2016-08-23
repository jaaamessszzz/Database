#!/usr/bin/python2.4
# encoding: utf-8
"""
primers.py

A script to automatically add primers to the Plasmids database from TSV files in user directories.

Created by Shane O'Connor 2016.
"""

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


def parse(dbi, primers_file, user_id, errors = []):

    # Import data
    df = pd.read_csv(primers_file, delimiter = '\t', encoding = 'utf-8')

    # Parse and add to model
    lcase_columns = [cn.lower() for cn in df.columns.values]
    for f in Primers._required_fields:
        if f not in df.columns.values:
            errors.append('Missing required field "{0}" in {1}.'.format(f, primers_file))
    for f in Primers._derived_or_generated_fields:
        if f in df.columns.values or f.lower() in lcase_columns:
            errors.append('Found derived field "{0}" in {1}. Since these values are derived during import into the database, storing them explicitly in the file could lead to inconsistencies.'.format(f, primers_file))

    if errors:
        raise Exception('Failed to parse the TSV file {0}'.format(primers_file))

    # Create empty dataframes used to store the successfully imported records and failed records respectively
    pass_df, fail_df = pd.DataFrame(columns = df.columns.values), pd.DataFrame(columns = df.columns.values)

    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    for index, row in df.iterrows():

        # Set up derived fields. creator_entry_number is created automatically by a database trigger upon import
        d = dict(
            creator = user_id,
            date = timestamp,
        )
        for header in Primers._required_fields + Primers._optional_fields:
            d[header] = row.get(header, None)

        # Use a new transaction per row (we could relax this)
        tsession = dbi.get_session(new_session = True)

        try:
            Primers.add(tsession, d)
            tsession.commit()
            tsession.close()
            pass_df = pass_df.append(row, ignore_index = True)
        except Exception, e:
            errors.append('Error: {0} on record {1}'.format(str(e), str(d)))
            fail_df = fail_df.append(row, ignore_index = True)
            tsession.rollback()
            tsession.close()

    dt = datetime.datetime.now()
    archive_file_path = os.path.join(os.path.split(primers_file)[0], '.' + os.path.splitext(os.path.split(primers_file)[1])[0] + '.{0}.archive.tsv'.format(dt.strftime('%Y%m%d_%H%M_%S_%f')))

    pass_df.to_csv(archive_file_path, sep = '\t', header = True, index = False)
    fail_df.to_csv(primers_file, sep = '\t', header = True, index = False)



def main():

    # Create up the database session
    dbi = DatabaseInterface(can_email = True)
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
                case_errors = []
                try:
                    parse(dbi, primers_file, user_id, case_errors)
                    if case_errors:
                        errors.append("Errors occurred processing '{0}':\n\t{1}".format(primers_file, '\n\t'.join(case_errors)))
                        colortext.warning(errors[-1])
                except Exception, e:
                    errors.append("Errors occurred processing '{0}': {1}\n\t{2}\n{3}".format(primers_file, str(e), '\n\t'.join(case_errors), traceback.format_exc()))
                    colortext.warning('Error: {0}\n{1}'.format(str(e), traceback.format_exc()))

    if errors:
        colortext.warning('\n\n'.join(errors))
        dbi.email_plasdmins('Automatic parsing of primers failed.', '<br/><br/>'.join([m.replace('\n', '<br/>') for m in errors]), '\n\n'.join(errors))



if __name__ == '__main__':
    main()





