#!/usr/bin/python2.4
# encoding: utf-8
"""
interface.py

Plasmids database interface code.

Created by Shane O'Connor 2016.
"""

import os
import json
import pprint
import traceback

import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText

from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import create_engine, and_

from klab import colortext
from klab.fs.fsio import read_file
from klab.db.sqlalchemy_interface import get_or_create_in_transaction, MySQLSchemaConverter

import model
from model import declarative_base, DBConstants


class DatabaseInterface(object):

    def __init__(self, echo_sql = False, can_email = True, config_file = None):

        # Read config
        if not config_file:
            config_file = os.path.expanduser('~/.my.cnf.json')
        try:
            config = json.loads(read_file(config_file))
        except:
            config = json.loads(read_file('/usr/local/turbogears/lablocal/.my.cnf.json'))
        self.config = config
        self.echo_sql = echo_sql
        self.can_email = can_email

        # Set up SQLAlchemy connections
        try:
            self.connect_string = config['plasmids']['connect_string']
            self.connect_string_utf = config['plasmids']['connect_string_utf']
        except:
            pass
        self.engine, self.session = None, None
        self.engine_utf, self.session_utf = None, None
        self.get_engine(utf = False)
        self.get_engine(utf = True)
        self.get_session(utf = False)
        self.get_session(utf = True)

        #hostname = sys_settings.database.hostname, rosetta_scripts_path = Nonerosetta_database_path = None, cache_dir = sys_settings.cache.cache_dir, echo_sql = False, port = sys_settings.database.port, file_content_buffer_size = None):
        #, passwd, connect_string, connect_string_utf, username = sys_settings.database.username,


    #############################
    #                           #
    #  SQLAlchemy Engine setup  #
    #                           #
    #############################



    def get_engine(self, utf = False):
        if utf:
            if not self.engine_utf:
                self.engine_utf = create_engine(self.connect_string_utf, echo = self.echo_sql)
            return self.engine_utf
        else:
            if not self.engine:
                self.engine = create_engine(self.connect_string, echo = self.echo_sql)
            return self.engine


    def get_connection(self, utf = False):
        '''Used for raw MySQL queries
            e.g.
                connection = dbi.get_connection();
                for r in connection.execute("SELECT * FROM Primers"):
                    # r is a ResultProxy object
                    print(r.sequence)

       '''
        engine = self.get_engine(utf = utf)
        return engine.connect()


    def get_session(self, new_session = False, autoflush = True, autocommit = False, utf = False):
        engine = self.get_engine(utf = utf)
        if new_session or ((not(utf) and not(self.session)) or (utf and not(self.session_utf))):
            database_session_maker = sessionmaker(autoflush = autoflush, autocommit = autocommit)
            s = scoped_session(database_session_maker)
            DeclarativeBasePlasmids = declarative_base()
            metadata_plasmids = DeclarativeBasePlasmids.metadata
            s.configure(bind=engine)
            metadata_plasmids.bind = engine
            if new_session:
                return s
            else:
                if utf:
                    self.session_utf = s
                else:
                    self.session = s
        if utf:
            return self.session_utf
        else:
            return self.session


    def renew(self, utf = False):
        self.session = self.get_session(new_session = True, utf = utf)


    #############################
    #                           #
    #  Data retrieval           #
    #                           #
    #############################


    def get_admin_contacts(self):
        return [r.Value for r in self.get_session().query(DBConstants).filter(DBConstants.ParameterType == u'admin_contact').order_by(DBConstants.Parameter)]


    #############################
    #                           #
    #  Admin                    #
    #                           #
    #############################


    def email_plasdmins(self, subject, htmltext, plaintext):
        if self.can_email:
            recipients = ['shane.oconnor@ucsf.edu']
            if recipients:
                recipients = sorted(set(recipients))
                msg = MIMEMultipart('alternative')
                msg['Subject'] = subject
                msg['From'] = 'support@kortemmelab.ucsf'
                msg['To'] = ';'.join(recipients)
                msg['Reply-To'] = 'support@kortemmelab.ucsf'

                part1 = MIMEText(plaintext, 'plain', "utf-8")
                part2 = MIMEText(htmltext, 'html', "utf-8")
                msg.attach(part1)
                msg.attach(part2)

                s = smtplib.SMTP()
                s.connect()
                s.sendmail(msg['From'], recipients, msg.as_string())
                s.close()

# Example transaction code for cassette plasmid insertion
#
#tsession = dbi.get_session(new_session = True)
#try:
#    ...
#  <either create a main record using the klab package>
#    cassette_plasmid = get_or_create_in_transaction(tsession, model, values, missing_columns = [], variable_columns = [], updatable_columns = [], only_use_supplied_columns = False)
#  <or manually>
#    cassette_plasmid = Primers(**dict_with_field_values)
#    tsession.add(cassette_plasmid)
#    tsession.flush()
#    ...
#    <create all records needed for a cassette plasmid>
#    ...
#    tsession.commit()
#    tsession.close()
#    return plasmid_record
#except:
#    colortext.error('Failure.')
#    tsession.rollback()
#    tsession.close()
#    raise


def test():
    # Test
    dbi = DatabaseInterface()
    q = dbi.get_session().query(model.Primers)
    for p in q:
        print(p)
        print('')
    connection = dbi.get_connection();
    for r in connection.execute("SELECT * FROM Primers"):
        print('{0} {1} {2}'.format(r.creator, r.creator_entry_number, r.sequence))


if __name__ == '__main__':
    dbi = DatabaseInterface()
    config = dbi.config['plasmids']
    pprint.pprint(config)
    sc = MySQLSchemaConverter(config['user'], config['host'], config['database'], config['password'], config['port'], config['socket'])
    sc.get_sqlalchemy_schema()
