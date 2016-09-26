import pprint
import re
import traceback

import numpy

from sqlalchemy import Table, ForeignKey, Column
from sqlalchemy.dialects.mysql import DOUBLE, TINYINT, LONGBLOB
from sqlalchemy.types import DateTime, Enum, Float, TIMESTAMP, Integer, Text, Unicode, String
from sqlalchemy.orm import relationship

# declarative_base instance used for TurboGears integration
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import inspect as sqlalchemy_inspect
from sqlalchemy import and_, or_, func
DeclarativeBasePlasmid = declarative_base()

from klab import colortext
from klab.bio.basics import translate_codons


# todo: this should be moved into the klab repository e.g. klab/bio/dna.py
class NotDNAException(Exception): pass
def reverse_complement(input_sequence):
    base_pair = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
        }
    try:
        return ''.join(reversed([base_pair[base.upper()] for base in input_sequence]))
    except:
        raise NotDNAException('The sequence is not a DNA (ACGT) sequence.')


def row_to_dict(r, keep_relationships = False):
    '''Converts an SQLAlchemy record to a Python dict. We assume that _sa_instance_state exists and is the only value we do not care about.
       If DeclarativeBase is passed then all DeclarativeBase objects (e.g. those created by relationships) are also removed.
    '''
    d = {}
    if not keep_relationships:
        # only returns the table columns
        t = r.__table__
        for c in [c.name for c in list(sqlalchemy_inspect(t).columns)]:
            d[c] = getattr(r, c)
        return d
    else:
        # keeps all objects including those of type DeclarativeBase or InstrumentedList and the _sa_instance_state object
        return copy.deepcopy(r.__dict__)


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


    def get_id(self):
        '''Return a canonically formatted string to identify the primer in stores e.g. "oJL0023".'''
        return 'o{0}{1:04d}'.format(self.creator, self.creator_entry_number)


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

            return db_record_object
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
    plasmid_type = Column(Enum('part', 'cassette', 'multicassette', 'other'), nullable=False)
    location = Column(Unicode(250, collation="utf8_bin"), nullable=False)
    description = Column(Text(), nullable=False)
    sequence = Column(Text(), nullable=False)
    status = Column(Enum('designed', 'verified', 'abandoned'), nullable=False)
    date = Column(TIMESTAMP, nullable=False)


    def get_id(self):
        '''Return a canonically formatted string to identify the plasmid in stores e.g. "pJL0023".'''
        return 'p{0}{1:04d}'.format(self.creator, self.creator_entry_number)


    def get_details(self, tsession):
        '''.'''
        d = row_to_dict(self)

        # Store the generated ID
        d['id'] = self.get_id()

        # Add a list to inform the user about errors or inconsistencies in the records
        d['errors'] = []

        # Rename plasmid_name to match the Primers table
        d['secondary_id'] = d['plasmid_name']
        del d['plasmid_name']

        # Retrieve the list of associated files
        d['files'] = []
        plasmid_files = tsession.query(Plasmid_File).filter(and_(Plasmid_File.creator == self.creator, Plasmid_File.creator_entry_number == self.creator_entry_number)).order_by(Plasmid_File.file_name)
        if plasmid_files.count() > 0:
            for pf in plasmid_files:
                d['files'].append(pf.get_details(tsession))

        # Retrieve the list of associated features
        d['features'] = []
        plasmid_features = tsession.query(Plasmid_Feature).filter(and_(Plasmid_Feature.creator == self.creator, Plasmid_Feature.creator_entry_number == self.creator_entry_number)).order_by(Plasmid_Feature.feature_name)
        if plasmid_features.count() > 0:
            for pf in plasmid_features:
                d['features'].append(pf.get_details(tsession))


        # todo: add all relevant information to part index list (only name and indicies so far)
        part_indices_list = None
        if self.plasmid_type.lower() == 'cassette':

            part_indices_list = []
            # List of [part_name, (starting_index, ending_index)]

            cassette_parts = tsession.query(Cassette_Assembly, Plasmid) \
                .filter(and_(self.creator == Cassette_Assembly.Cassette_creator,
                             self.creator_entry_number == Cassette_Assembly.Cassette_creator_entry_number,
                             Plasmid.creator == Cassette_Assembly.Part_creator,
                             Plasmid.creator_entry_number == Cassette_Assembly.Part_creator_entry_number))

            print('*!' * 100)
            print(cassette_parts)
            for cassette_assembly, plasmid in cassette_parts:
                # Do the restriction digest with BsaI and get indices for start/end of part minus overhangs
                print('***')
                print(row_to_dict(cassette_assembly))
                part_sequence = self.fetch_cassette_parts(tsession, [(plasmid.creator, plasmid.creator_entry_number)])
                if len(part_sequence) == 1:
                    part_sequence = part_sequence[0][4:-4]
                    for instance in re.finditer(part_sequence.upper().strip(), self.sequence):
                        part_indices_list.append([plasmid.description, cassette_assembly.Part_number, (instance.start(), instance.end())])

        d['part_indices'] = part_indices_list
        print("d['part_indices']", d['part_indices'])

        # Retrieve the list of plasmid type-specific details
        d['details'] = None
        d['printed_plasmid_type'] = d['plasmid_type'].title() # used by the website for longer descriptions e.g. "Part 3, 4"
        part_numbers = None
        if self.plasmid_type == 'part':
            resistance, part_number = None, None
            try:
                resistance = tsession.query(Part_Plasmid).filter(and_(Part_Plasmid.creator == self.creator, Part_Plasmid.creator_entry_number == self.creator_entry_number)).one().resistance
            except:
                d['errors'].append('This part plasmid is missing a Part_Plasmid record so the resistance is unknown.')
            try:
                plasmid_parts = tsession.query(Part_Plasmid_Part).filter(and_(Part_Plasmid_Part.creator == self.creator, Part_Plasmid_Part.creator_entry_number == self.creator_entry_number))
                part_numbers = sorted([r.part_number for r in plasmid_parts])
            except Exception, e:
                colortext.warning(str(e))
                colortext.warning(traceback.format_exc())
                d['errors'].append('This part plasmid is missing a Part_Plasmid_Part record so the part number is unknown.')

            d['printed_plasmid_type'] = '{0} {1}'.format(self.plasmid_type.title(), ', '.join(part_numbers))
            d['details'] = dict(
                resistance = resistance,
                part_numbers = part_numbers,
            )
        elif self.plasmid_type == 'cassette':
            assert ('part_plasmids' not in d)
            d['part_plasmids'] = []
            try:
                parts_read = {}
                part_plasmids = tsession.query(Cassette_Assembly).filter(and_(Cassette_Assembly.Cassette_creator == self.creator, Cassette_Assembly.Cassette_creator_entry_number == self.creator_entry_number)).order_by(Cassette_Assembly.Part_number)
                for pp in part_plasmids:
                    pp_id = (pp.Part_creator, pp.Part_creator_entry_number)
                    if pp_id not in parts_read:
                        parts_read[pp_id] = True
                        assert(pp.Part_number not in d['part_plasmids'])
                        part_plasmid_record = tsession.query(Plasmid).filter(and_(Plasmid.creator == pp.Part_creator, Plasmid.creator_entry_number == pp.Part_creator_entry_number)).one()
                        assert(part_plasmid_record.plasmid_type == 'part')
                        d['part_plasmids'].append(part_plasmid_record.get_details(tsession))
            except Exception, e:
                colortext.warning(str(e))
                colortext.warning(traceback.format_exc())
                d['errors'].append('An error occurred retrieving data for this cassette assembly.')

        d['details'] = d['details'] or {}

        return d



    def fetch_cassette_parts(self, tsession, input_sequences, table_info=None):
        # Get sequences and part numbers for listed part plasmids

        part_IDs = [and_(Plasmid.creator == part_creator, Plasmid.creator_entry_number == part_creator_entry_number) for (part_creator, part_creator_entry_number) in input_sequences]
        part_plasmids_query = tsession.query(Plasmid, Part_Plasmid, Part_Plasmid_Part, Part_Type) \
            .filter(and_(Part_Plasmid.creator == Plasmid.creator,
                         Part_Plasmid.creator_entry_number == Plasmid.creator_entry_number,
                         Part_Plasmid.creator == Part_Plasmid_Part.creator,
                         Part_Plasmid.creator_entry_number == Part_Plasmid_Part.creator_entry_number,
                         Part_Plasmid_Part.part_number == Part_Type.part_number,
                         or_(*part_IDs)
                         )
                    )

        restriction_enzyme = 'BsaI'
        site_F = 'GGTCTC'
        site_R = 'GAGACC'

        already_fetched_cassettes = []
        part_sequences = []

        for plasmid, part_plasmid, part_plasmid_part, part_type in part_plasmids_query:

            print "%s\t%s\t%s\t%s" % (
            plasmid.plasmid_name, plasmid.creator, plasmid.creator_entry_number, part_plasmid_part.part_number)

            sequence_upper = plasmid.sequence.upper()

            if sequence_upper.count(site_F) != 1:
                raise Plasmid_Exception('There is more than one forward %s site in %s! %s' % (restriction_enzyme, plasmid.plasmid_name or plasmid.get_id(), plasmid.sequence.count(site_F)))
            if sequence_upper.count(site_R) != 1:
                raise Plasmid_Exception('There is more than one reverse %s site in %s!' % (restriction_enzyme, plasmid.plasmid_name or plasmid.get_id()))

            left_overhang = part_type.overhang_5
            right_overhang = part_type.overhang_3

            # check if a part plasmid has multiple part plasmid parts - will select the correct 5' and 3' overhangs
            if part_plasmids_query.filter( Plasmid.creator == plasmid.creator).filter(Plasmid.creator_entry_number == plasmid.creator_entry_number).count() > 1:
                left_overhangs = set()
                right_overhangs = set()

                for p, pp, ppp, pt in part_plasmids_query.filter(and_(Plasmid.creator == plasmid.creator, Plasmid.creator_entry_number == plasmid.creator_entry_number)):
                    left_overhangs.add(pt.overhang_5)
                    right_overhangs.add(pt.overhang_3)

                leftoverhangs = left_overhangs.symmetric_difference(right_overhangs) # Leftover overhangs? Get it?

                for p, pp, ppp, pt in part_plasmids_query.filter(and_(Plasmid.creator == plasmid.creator, Plasmid.creator_entry_number == plasmid.creator_entry_number)):
                    if pt.overhang_5 in leftoverhangs:
                        left_overhang = pt.overhang_5
                    if pt.overhang_3 in leftoverhangs:
                        right_overhang = pt.overhang_3

            site_F_position = sequence_upper.find(left_overhang, sequence_upper.find(site_F))
            site_R_position = sequence_upper.find(right_overhang, sequence_upper.find(site_R) - 10) + 4

            # Circular permutation so that forward cut site is always upstream of reverse cut site in linear sequence
            if site_F_position > site_R_position:
                sequence_upper = sequence_upper[site_R_position + 8:] + sequence_upper[:site_R_position + 8]

            if table_info:
                if (plasmid.creator, plasmid.creator_entry_number) not in already_fetched_cassettes:
                    table_info['Part list'].append(sequence_upper[
                                                   sequence_upper.find(left_overhang, sequence_upper.find(site_F)):
                                                   sequence_upper.find(right_overhang,
                                                                       sequence_upper.find(site_R) - 10) + 4
                                                   ])
                    table_info['Description'].append(plasmid.description)
                    table_info['Part list ID'].append(
                        (part_plasmid_part.part_number, plasmid.creator, plasmid.creator_entry_number))
                    if part_plasmid_part.part_number == '1' or part_plasmid_part.part_number == '5':
                        table_info['Connectors'][part_plasmid_part.part_number] = [plasmid.creator,
                                                                                   plasmid.creator_entry_number]
                    already_fetched_cassettes.append((plasmid.creator, plasmid.creator_entry_number))

            part_sequences.append(sequence_upper[sequence_upper.find(left_overhang, sequence_upper.find(site_F)):sequence_upper.find(right_overhang,sequence_upper.find(site_R) - 10) + 4])

        if table_info:
            return table_info
        else:
            return part_sequences



    @staticmethod
    def add(tsession, input_dict, silent=True):

        try:
            for f in Plasmid._optional_fields:
                if f not in input_dict:
                    input_dict[f] = None

            db_record_object = Plasmid(**input_dict)

            tsession.add(db_record_object)
            tsession.flush()

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            # If the user did not specify a plasmid name, use the auto-generated name
            if not db_record_object.plasmid_name:
                db_record_object.plasmid_name = db_record_object.get_id()
            tsession.flush()

            return db_record_object
        except:
            raise


    def __repr__(self):
        return 'Internal numbering: {0} {1}\n' \
            'plasmid_name: {2}'.format(self.creator, self.creator_entry_number, self.plasmid_name)

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

            return db_record_object
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

            return db_record_object
        except:
            raise


class Feature(DeclarativeBasePlasmid):
    __tablename__ = 'Feature'

    Feature_name = Column(Unicode(200, collation="utf8_bin"), nullable=False, primary_key=True)
    MD5_hash = Column(String(64), nullable=True)
    Feature_type = Column(Unicode(200, collation="utf8_bin"), ForeignKey('Feature_Type.Feature_type'), nullable=False)
    Feature_sequence = Column(Text(collation="utf8_bin"), nullable=True)
    description = Column(Text(collation="utf8_bin"))

    ftype = relationship('Feature_Type', primaryjoin='Feature_Type.Feature_type == Feature.Feature_type')


    def get_details(self, tsession):
        color = tsession.query(Feature_Type).filter(Feature_Type.Feature_type == self.Feature_type).one().color

        protein_sequence_translation = None
        try:
            if len(self.Feature_sequence) % 3 == 0:
                protein_sequence_translation = translate_codons(self.Feature_sequence)
        except: pass

        return dict(
                name = self.Feature_name,
                hash = self.MD5_hash,
                type = self.Feature_type,
                sequence = self.Feature_sequence,
                protein_sequence = protein_sequence_translation,
                description = self.description,
                color = color,
        )


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


    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Multicassette_Assembly(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()

            return db_record_object
        except:
            raise

class Multicassette_Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Multicassette_Plasmid'

    creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    creator = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    Backbone_creator = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    Backbone_creator_entry_number = Column(Integer, nullable=False, primary_key=True)

    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Multicassette_Plasmid(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()

            return db_record_object
        except:
            raise

class Origin(DeclarativeBasePlasmid):
    __tablename__ = 'Origin'

    Origin = Column(Unicode(100, collation="utf8_bin"), nullable=False, primary_key=True)


class Part_Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Part_Plasmid'

    _required_fields = ['creator', 'creator_entry_number', 'resistance']

    creator_entry_number = Column(Integer, ForeignKey('Plasmid.creator_entry_number'), nullable=False, primary_key=True)
    creator = Column(Unicode(5), ForeignKey('Plasmid.creator'), nullable=False, primary_key=True)
    resistance = Column(Unicode(100, collation="utf8_bin"), nullable=False)

    plasmid = relationship('Plasmid', primaryjoin='and_(Plasmid.creator == Part_Plasmid.creator, Plasmid.creator_entry_number == Part_Plasmid.creator_entry_number)')
    parts = relationship('Part_Plasmid_Part', primaryjoin='and_(Part_Plasmid.creator == Part_Plasmid_Part.creator, Part_Plasmid.creator_entry_number == Part_Plasmid_Part.creator_entry_number)')


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

            return db_record_object
        except:
            raise


    def to_dict(self):
        '''This function is used by the web interface.'''
        d = row_to_dict(self)
        d.update(row_to_dict(self.plasmid))
        d['parts'] = sorted([p.part_number for p in self.parts])
        d['__id__'] = self.plasmid.get_id() # this is easier to work with for the website hashtables
        return d


    def __repr__(self):
        return 'Internal numbering: {0} {1}'.format(self.creator, self.creator_entry_number)


class Part_Plasmid_Part(DeclarativeBasePlasmid):
    __tablename__ = 'Part_Plasmid_Part'

    _required_fields = ['creator', 'creator_entry_number', 'part_number']

    creator_entry_number = Column(Integer, ForeignKey('Part_Plasmid.creator_entry_number'), nullable=False, primary_key=True)
    creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Part_Plasmid.creator'), nullable=False, primary_key=True)
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

            return db_record_object
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
    creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Plasmid.creator'), nullable=False)
    creator_entry_number = Column(Integer, ForeignKey('Plasmid.creator_entry_number'), nullable = False)
    feature_name = Column(Unicode(200, collation="utf8_bin"), ForeignKey('Feature.Feature_name'), nullable=False)

    feature = relationship('Feature', viewonly = True, primaryjoin = "Plasmid_Feature.feature_name==Feature.Feature_name")
    plasmid = relationship('Plasmid', viewonly = True, primaryjoin = 'and_(Plasmid.creator == Plasmid_Feature.creator, Plasmid.creator_entry_number == Plasmid_Feature.creator_entry_number)')


    def get_details(self, tsession):
        '''Returns details about a feature.'''
        feature_details = self.feature.get_details(tsession)
        plasmid_sequence = self.plasmid.sequence
        feature_sequence = self.feature.Feature_sequence.upper().strip()
        feature_type = self.feature.ftype
        indices = []

        for instance in re.finditer(feature_sequence, plasmid_sequence):
            if self.feature.Feature_type == 'Restrxn Type II':
                indices.append(dict(
                    name = self.feature_name + ' F',
                    start = instance.start(),
                    end = instance.end(),
                ))
                indices.append(dict(
                    name = self.feature_name + ' Overhang',
                    start = instance.start() + 7,
                    end = instance.start() + 11,
                ))
            else:
                indices.append(dict(
                        name = self.feature_name,
                        start = instance.start(),
                        end = instance.end(),
                ))

        for instance in re.finditer(reverse_complement(feature_sequence), plasmid_sequence):
            if self.feature.Feature_type == 'Restrxn Type II':
                indices.append(dict(
                        name = self.feature_name + ' R',
                        start = instance.start(),
                        end = instance.end(),
                ))
                indices.append(dict(
                        name = self.feature_name + ' Overhang',
                        start = instance.start() - 5,
                        end = instance.start() - 1,
                ))
            else:
                indices.append(dict(
                        name = self.feature_name,
                        start = instance.start(),
                        end = instance.end(),
                ))
        feature_details['indices'] = indices
        return feature_details


    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Plasmid_Feature(**input_dict)
            if not silent:
                colortext.pcyan('Adding this record:')
                print('')
            feature_check = tsession.query(Plasmid_Feature).filter(and_(Plasmid_Feature.creator == db_record_object.creator, Plasmid_Feature.creator_entry_number == db_record_object.creator_entry_number, Plasmid_Feature.feature_name == db_record_object.feature_name))
            if feature_check.count() == 0:
                tsession.add(db_record_object)
                tsession.flush()
                return db_record_object
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
    color = Column(String(6), nullable = False)


class Plasmid_File(DeclarativeBasePlasmid):
    __tablename__ = 'Plasmid_File'

    _required_fields = ['creator', 'creator_entry_number', 'file_name', 'file_type', 'Description', 'File']

    ID = Column(Integer, nullable=False, primary_key=True)
    creator_entry_number = Column(Integer, ForeignKey('Part_Plasmid_Part.creator_entry_number'), nullable=False)
    creator = Column(Unicode(5), ForeignKey('Part_Plasmid_Part.creator'), nullable=False)
    file_name = Column(Unicode(100, collation="utf8_bin"), nullable=False)
    file_type = Column(Unicode(100), nullable=False)
    Description = Column(Text(), nullable=False)
    File = Column(LONGBLOB(), nullable=False)


    def get_details(self, tsession):
        '''Returns meta-data about a file.'''

        # todo: For readability, I should really be defining relationships above and using them here
        ft = tsession.query(File_Type).filter(File_Type.file_type == self.file_type).one()
        file_extension = ft.file_extension
        file_color = ft.color

        return dict(
            ID = self.ID,
            creator_entry_number = self.creator_entry_number,
            creator = self.creator,
            file_name = self.file_name,
            file_type = self.file_type,
            file_extension = file_extension,
            Description = self.Description,
            length_in_bytes = len(self.File),
            color = file_color,
        )


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

            return db_record_object
        except:
            raise


    def __repr__(self):
        return '{0}, type {1}, created by {2}: {3}'.format(self.file_name, self.file_type, self.creator, self.Description)

class CDS_Mutant(DeclarativeBasePlasmid):
    __tablename__ = 'CDS_Mutant'

    ID = Column(Integer, nullable=False, primary_key=True)
    Plasmid_Feature_ID = Column(Integer, ForeignKey('Plasmid_Feature.ID'), nullable=False)
    Mutant_ID = Column(Integer, nullable=False)
    creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Users.ID'), nullable=False)
    Description = Column(Text(), nullable=True)
    date = Column(TIMESTAMP, nullable=False)

    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = CDS_Mutant(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()

            return db_record_object
        except:
            raise

class CDS_Mutant_Constituent(DeclarativeBasePlasmid):
    __tablename__ = 'CDS_Mutant_Constituent'

    ID = Column(Integer, ForeignKey('CDS_Mutant.ID'), nullable=False, primary_key=True)
    mutation = Column(Unicode(16, collation="utf8_bin"), nullable=False)
    position = Column(Integer, nullable=False)
    wt_AA = Column(Unicode(1, collation="utf8_bin"), nullable=False)
    mut_AA = Column(Unicode(1, collation="utf8_bin"), nullable=False)
    Description = Column(Text(), nullable=True)

    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = CDS_Mutant_Constituent(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()

            return db_record_object
        except:
            raise