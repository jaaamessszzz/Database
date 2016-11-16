import pprint
import re
import traceback

import numpy

from sqlalchemy import Table, ForeignKey, Column, UniqueConstraint, FetchedValue
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

from procedures import call_procedure


# todo: this should be moved into the klab repository e.g. klab/bio/dna.py
class NotDNAException(Exception): pass



class DuplicatePublicationException(Exception): pass

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
        raise NotDNAException('The sequence is not a DNA (AGCT) sequence.')


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
    secondary_id = Column(Unicode(100, collation = "utf8_bin"), nullable = True)
    sequence = Column(Unicode(256, collation="utf8_bin"), nullable = False)
    direction = Column(Enum('F','R'), nullable=False)
    description = Column(Text(collation="utf8_bin"), nullable=False)
    TM = Column(Float, nullable=False)
    date = Column(TIMESTAMP, nullable=False)
    GC = Column(Float, nullable=True)
    length = Column(Float, nullable=True)
    template_id = Column(Unicode(200, collation="utf8_bin"), nullable = True)
    template_description = Column(Text(collation="utf8_bin"), nullable = True)

    UniqueConstraint('sequence', 'sequence', name = 'secondary_id_2')
    UniqueConstraint('sequence', 'direction', name = 'secondary_id_3')


    def get_id(self):
        '''Return a canonically formatted string to identify the primer in stores e.g. "oJL0023".'''
        return 'o{0}{1:04d}'.format(self.creator, self.creator_entry_number)


    def to_dict(self):
        '''This function is used by the web interface.'''
        d = row_to_dict(self)
        d['__id__'] = self.get_id() # this is easier to work with for the website hashtables
        return d


    @staticmethod
    def add(tsession, d, silent = True):

        try:
            sequence = d['sequence']

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

            assert(len(sequence) < 256)

            # Calculate GC content
            gc = Primers.calculate_gc(sequence)
            if 'GC' in d and d['GC'] != None:
                assert(abs(d['GC'] - gc) < 0.01)
            d['GC'] = gc

            # Calculate length
            primer_length = len(sequence)
            if 'length' in d and d['length'] != None:
                assert(d['length'] == primer_length)
            d['length'] = primer_length

            db_record_object = Primers(**d)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()

            # Note: Because we use a MySQL trigger, SQLAlchemy does not know to update creator_entry_number so db_record_object.creator_entry_number will return zero
            #       I looked into this for an hour and tried using FetchedValue but gave up on it. Other solutions: i) use PostgreSQL; ii) define the trigger using DDL in SQLAlchemy.

            if not silent:
                colortext.pcyan('Added this record:')
                print(db_record_object)
                print('')

            return db_record_object
        except:
            raise


    @staticmethod
    def calculate_gc(sequence):
        gc = 0
        for c in sequence.lower():
            if c == 'g' or c == 'c':
                gc += 1
            else:
                assert(c == 'a' or c == 't')
        return float(gc) / float(len(sequence))


    def add_primer_plasmid(self, tsession, associated_plasmid, silent = True):
        '''
        Add a PlasmidPrimer record for this primer
        This was meant to be a staticmethod which added both the Primers record and the PlasmidPrimer record.
        However, I hit a snag. Because we use a MySQL trigger, SQLAlchemy does not know to update creator_entry_number
        so db_record_object.creator_entry_number will return zero. There are potential fixes but I looked into this
        for an hour and gave up on it.

        Potential fixes: i) use PostgreSQL; ii) define the trigger using DDL in SQLAlchemy.

        :param tsession:
        :param associated_plasmid:
        :param silent:
        :return:
        '''
        try:
            # Add PlasmidPrimer record
            db_record_object = PlasmidPrimer(**dict(
                plasmid_creator = associated_plasmid.creator,
                plasmid_creator_entry_number = associated_plasmid.creator_entry_number,
                primer_creator = self.creator,
                primer_creator_entry_number = self.creator_entry_number,
            ))

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
    suitable_for_parts = Column(TINYINT(1), default = 0, nullable = False)
    suitable_for_others = Column(TINYINT(1), default = 1, nullable = False)
    suitable_for_templates = Column(TINYINT(1), default = 1, nullable = False)


class Users(DeclarativeBasePlasmid):
    __tablename__ = 'Users'

    ID = Column(Unicode(5, collation="utf8_bin"), nullable=False, primary_key=True)
    github_username = Column(Unicode(26, collation="utf8_bin"), nullable=False)
    lab_username = Column(Unicode(50, collation="utf8_bin"), nullable = False)
    color = Column(Unicode(6, collation="utf8_bin"), nullable = False)


class Publication_Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Publication_Plasmid'

    ID = Column(Integer, nullable = False, primary_key = True, autoincrement = True)
    PublicationID = Column(Integer, ForeignKey('Publication.ID'), nullable = False)
    creator = Column(Unicode(5), ForeignKey('Plasmid.creator'), nullable = False)
    creator_entry_number = Column(Integer, ForeignKey('Plasmid.creator_entry_number'), nullable = False)
    Description = Column(Text, nullable = True)
    Deprecated = Column(TINYINT(1), default = 0, nullable = False)


class PlasmidPrimer(DeclarativeBasePlasmid):
    __tablename__ = 'PlasmidPrimer'

    plasmid_creator = Column(Unicode(5), ForeignKey('Plasmid.creator'), nullable = False, primary_key = True)
    plasmid_creator_entry_number = Column(Integer, ForeignKey('Plasmid.creator_entry_number'), nullable = False, primary_key = True)
    primer_creator = Column(Unicode(5), ForeignKey('Primers.creator'), nullable = False, primary_key = True)
    primer_creator_entry_number = Column(Integer, ForeignKey('Primers.creator_entry_number'), nullable = False, primary_key = True)

    plasmid = relationship('Plasmid', primaryjoin = 'and_(Plasmid.creator == PlasmidPrimer.plasmid_creator, Plasmid.creator_entry_number == PlasmidPrimer.plasmid_creator_entry_number)')
    primer = relationship('Primers', primaryjoin = 'and_(Primers.creator == PlasmidPrimer.primer_creator, Primers.creator_entry_number == PlasmidPrimer.primer_creator_entry_number)')


    def __repr__(self):
        return 'Plasmid: {0}, Primer: {1}'.format(self.plasmid, self.primer)


class Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Plasmid'

    _required_fields = ['sequence', 'description', 'location', 'plasmid type', 'status']
    _optional_fields = ['plasmid_name']
    _derived_or_generated_fields = ['creator', 'creator_entry_number', 'date']

    creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Users.ID'), nullable=False, primary_key=True)
    creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    plasmid_name = Column(Unicode(100, collation="utf8_bin"), nullable=True, unique=True)
    plasmid_type = Column(Enum('part', 'cassette', 'multicassette', 'other', 'template', 'design'), nullable=False)
    location = Column(Unicode(250, collation="utf8_bin"), nullable=False)
    description = Column(Text(collation="utf8_bin"), nullable=False)
    sequence = Column(Text(), nullable=False)
    status = Column(Enum('designed', 'verified', 'abandoned'), nullable=False)
    date = Column(TIMESTAMP, nullable=False)

    publications = relationship("Publication",
                                secondary = Publication_Plasmid.__table__,
                                primaryjoin = "and_(Plasmid.creator == Publication_Plasmid.creator, Plasmid.creator_entry_number == Publication_Plasmid.creator_entry_number, Publication_Plasmid.Deprecated == 0)",
                                secondaryjoin = "Publication_Plasmid.PublicationID == Publication.ID",
                                viewonly = True)

    primers = relationship("Primers",
                                secondary = PlasmidPrimer.__table__,
                                primaryjoin = "and_(Plasmid.creator == PlasmidPrimer.plasmid_creator, Plasmid.creator_entry_number == PlasmidPrimer.plasmid_creator_entry_number)",
                                secondaryjoin = "and_(PlasmidPrimer.primer_creator == Primers.creator, PlasmidPrimer.primer_creator_entry_number == Primers.creator_entry_number)",
                                viewonly = True)


    @staticmethod
    def get_by_id(tsession, creator, creator_entry_number):
        '''Return a canonically formatted string to identify the plasmid in stores e.g. "pJL0023".'''
        return tsession.query(Plasmid).filter(and_(Plasmid.creator == creator, Plasmid.creator_entry_number == creator_entry_number)).one()


    def get_id(self):
        '''Return a canonically formatted string to identify the plasmid in stores e.g. "pJL0023".'''
        return u'p{0}{1:04d}'.format(self.creator, self.creator_entry_number)


    @classmethod
    def get_plasmid_helper(cls, tsession, cache, creator, creator_entry_number):
        if cache:
            return cache.get_plasmid(creator, creator_entry_number)
        else:
            return tsession.query(cls).filter(and_(cls.creator == creator, cls.creator_entry_number == creator_entry_number)).one()


    ###############################
    # Generic plasmid properties
    ###############################


    def _get_details_basic(self):
        d = row_to_dict(self)

        # Store the generated ID
        d['id'] = self.get_id()

        # Add a list to inform the user about errors or inconsistencies in the records
        d['errors'] = []

        # Used for part type-specific properties
        d['details'] = {}

        # Used by the website for longer descriptions e.g. "Part 3, 4"
        d['printed_plasmid_type'] = d['plasmid_type'].title()

        # Rename plasmid_name to match the Primers table
        d['secondary_id'] = d['plasmid_name']
        del d['plasmid_name']
        return d


    ###############################
    # Detailed plasmid properties
    ###############################


    def _get_details_files(self, tsession, engine, d, cache = None, fresh = False):
        d = d or self._get_details_basic()
        plasmid_files = tsession.query(Plasmid_File).filter(and_(Plasmid_File.creator == self.creator, Plasmid_File.creator_entry_number == self.creator_entry_number)).order_by(Plasmid_File.file_name)
        if plasmid_files.count() > 0:
            for pf in plasmid_files:
                d['files'].append(pf.get_details(tsession, engine))


    def _get_details_publications(self, tsession, engine, d, cache = None, fresh = False):
        d = d or self._get_details_basic()
        for pub in self.publications:
            intersection_record = tsession.query(Publication_Plasmid).filter(and_(
                    Publication_Plasmid.creator == self.creator,
                    Publication_Plasmid.creator_entry_number == self.creator_entry_number,
                    Publication_Plasmid.PublicationID == pub.ID,
                    Publication_Plasmid.Deprecated == 0))
            if intersection_record:
                intersection_record = intersection_record.one()
                description = intersection_record.Description
                pubd = pub.get_details()
                pubd['description'] = description
                d['publications'].append(pubd)


    def _get_details_primers(self, tsession, engine, d, cache = None, fresh = False):
        d = d or self._get_details_basic()
        for prmr in self.primers:
            d['primers'].append(prmr.to_dict())


    def _get_details_features(self, tsession, engine, d, cache = None, fresh = False):
        d = d or self._get_details_basic()

        plasmid_features = tsession.query(Plasmid_Feature).filter(and_(Plasmid_Feature.creator == self.creator, Plasmid_Feature.creator_entry_number == self.creator_entry_number)).order_by(Plasmid_Feature.feature_name)
        if plasmid_features.count() > 0:
            for pf in plasmid_features:
                if cache:
                    d['features'].append(cache.get_feature(pf.feature_name).get_details(tsession, engine, cache = cache))
                else:
                    d['features'].append(pf.get_details(tsession, engine))


    def _get_details_part_indices_list(self, tsession, engine, d, cache = None, fresh = False, only_basic_details = False):
        '''Retrieve a List of dicts with [part_name, part_number, starting_index, ending_index].'''

        d = d or self._get_details_basic()
        cassette_parts = tsession.query(Cassette_Assembly, Plasmid) \
            .filter(and_(self.creator == Cassette_Assembly.Cassette_creator,
                         self.creator_entry_number == Cassette_Assembly.Cassette_creator_entry_number,
                         Plasmid.creator == Cassette_Assembly.Part_creator,
                         Plasmid.creator_entry_number == Cassette_Assembly.Part_creator_entry_number))

        # print('*!' * 100)
        # print(cassette_parts)
        for cassette_assembly, cassette_part_plasmid in cassette_parts:
            # Do the restriction digest with BsaI and get indices for start/end of part minus overhangs
            part_sequence, part_types = self.fetch_plasmid_parts(tsession, engine, cassette_part_plasmid.creator, cassette_part_plasmid.creator_entry_number, cache = cache, fresh = fresh)

            # todo: write an explanation here
            if len(part_sequence) == 1:
                part_sequence = part_sequence[0][4:-4]
                for instance in re.finditer(part_sequence.upper().strip(), self.sequence):
                    d['part_indices'].append(dict(
                        name = cassette_part_plasmid.description,
                        part_number = ', '.join(part_types),
                        start = instance.start(),
                        end = instance.end(),
                    ))


    ##################################
    # Plasmid type-specific properties
    ##################################


    def _get_cassette_plasmid_details(self, tsession, engine, d, cache = None, fresh = False, only_basic_details = False):

        d = d or self._get_details_basic()
        assert (self.plasmid_type == 'cassette')

        assert ('part_plasmids' not in d)
        d['part_plasmids'] = []
        if not only_basic_details:
            try:
                parts_read = {}
                part_plasmids = tsession.query(Cassette_Assembly).filter(
                    and_(Cassette_Assembly.Cassette_creator == self.creator,
                         Cassette_Assembly.Cassette_creator_entry_number == self.creator_entry_number)).order_by(
                    Cassette_Assembly.Part_number)
                for pp in part_plasmids:
                    pp_id = (pp.Part_creator, pp.Part_creator_entry_number)
                    if pp_id not in parts_read:
                        parts_read[pp_id] = True
                        assert (pp.Part_number not in d['part_plasmids'])

                        if cache:
                            d['part_plasmids'].append(cache.get_plasmid(pp.Part_creator, pp.Part_creator_entry_number).get_details(tsession, engine, cache = cache))
                            assert(d['part_plasmids'][-1]['plasmid_type'] == 'part')
                        else:
                            part_plasmid_record = tsession.query(Plasmid).filter(and_(Plasmid.creator == pp.Part_creator, Plasmid.creator_entry_number == pp.Part_creator_entry_number)).one()
                            assert (part_plasmid_record.plasmid_type == 'part')
                            d['part_plasmids'].append(part_plasmid_record.get_details(tsession, engine))




            except Exception, e:
                colortext.warning(str(e))
                colortext.warning(traceback.format_exc())
                d['errors'].append('An error occurred retrieving data for this cassette assembly.')


    def _get_other_plasmid_details(self, tsession, engine, d, cache = None, fresh = False, only_basic_details = False):

        d = d or self._get_details_basic()
        assert(self.plasmid_type == 'other')

        sub_record, resistance, vector = None, None, None
        try:
            sub_record = tsession.query(Other_Plasmid).filter(and_(Other_Plasmid.creator == self.creator, Other_Plasmid.creator_entry_number == self.creator_entry_number)).one()
        except:
            d['errors'].append('This plasmid is missing am Other_Plasmid record so the resistance is unknown.')
        if sub_record:
            # Storing the vector ID is optional
            resistance, vector = sub_record.resistance, sub_record.vector or None
        d['details'] = dict(resistance = resistance, vector = vector)


    def _get_design_plasmid_details(self, tsession, engine, d, cache = None, fresh = False, only_basic_details = False):

        d = d or self._get_details_basic()
        assert(self.plasmid_type == 'design')

        sub_record, resistance, vector, is_template = None, None, None, None
        try:
            sub_record = tsession.query(Design_Plasmid).filter(and_(Design_Plasmid.creator == self.creator, Design_Plasmid.creator_entry_number == self.creator_entry_number)).one()
        except:
            d['errors'].append('This plasmid is missing a Design_Plasmid record so the resistance, vector, and template status are unknown.')
        if sub_record:
            # Storing the vector ID is optional
            resistance, vector, is_template = sub_record.resistance, sub_record.vector or None, sub_record.is_template
        d['details'] = dict(resistance = resistance, vector = vector, is_template = is_template)


    def _get_part_plasmid_details(self, tsession, engine, d, cache = None, fresh = False, only_basic_details = False):

        d = d or self._get_details_basic()
        assert (self.plasmid_type == 'part')

        part_numbers = None
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


    ##################################
    # Main get_details function
    ##################################


    def get_details(self, tsession, engine, cache = None, fresh = False, only_basic_details = False):
        '''.'''

        if fresh or '__details__' not in self.__dict__:

            d = self._get_details_basic()

            # Fill in optional details
            detail_functions = dict(
                all = dict(
                    files = self._get_details_files,               # Retrieve the list of associated files
                    publications = self._get_details_publications, # Retrieve the list of associated publications
                    primers = self._get_details_primers,           # Retrieve the list of associated primers
                    features = self._get_details_features,         # Retrieve the list of associated features
                ),
                cassette = dict(
                    part_indices = self._get_details_part_indices_list,  # todo: add all relevant information to part index list (only name and indicies so far)
                    details = self._get_cassette_plasmid_details,
                ),
                other = dict(
                    details = self._get_other_plasmid_details,
                ),
                part = dict(
                    details = self._get_part_plasmid_details,
                ),
                design = dict(
                    details = self._get_design_plasmid_details,
                ),
            )

            # Fill in generic optional details
            for detail_type, f in sorted(detail_functions['all'].iteritems()):
                d[detail_type] = []
                if not only_basic_details:
                    f(tsession, engine, d, cache = cache, fresh = fresh)

            # Fill in plasmid type-specific details
            for detail_type, g in sorted(detail_functions.get(self.plasmid_type, {}).iteritems()):
                d[detail_type] = []
                if not only_basic_details:
                    g(tsession, engine, d, cache = cache, fresh = fresh, only_basic_details = only_basic_details)

            self.__dict__['__details__'] = d

        return self.__dict__['__details__']


    def fetch_plasmid_parts(self, tsession, engine, creator, creator_entry_number, cache = None, fresh = False):
        '''

        :param tsession:
        :param engine:
        :param input_plasmids ([creator, creator_entry_number]): The function accepts multiple input plasmids (Q: should it? maybe for gg assembly but not for get_details)
        :return:
        '''
        # Get sequences and part numbers for listed part plasmids

        restriction_enzyme = 'BsaI'
        site_F = 'GGTCTC'
        site_R = 'GAGACC'

        # todo: use caching here
        #if cache:
        #    d['features'].append(cache.get_feature(pf.feature_name).get_details(tsession, engine, cache = cache))
        #else:
        #    d['features'].append(pf.get_details(tsession, engine))

        # Retrieve the part plasmid record
        part_plasmid = self.get_plasmid_helper(tsession, cache, creator, creator_entry_number)
        assert(part_plasmid.plasmid_type == 'part')

        # Return the list of parts in the plasmid
        plasmid_parts = call_procedure(engine, 'getPartPlasmidParts', [creator, creator_entry_number], as_dict = True)

        assert(len(plasmid_parts) > 0)
        sequence = set([p['Plasmid_sequence_UC'] for p in plasmid_parts]) # these should all be the same - we could assert that
        assert(len(sequence) == 1)
        sequence = sequence.pop()
        part_sequences = [sequence] # todo: change to only expect one sequence

        part_types = [p['Part_Type_part_number'] for p in plasmid_parts]

        if sequence.count(site_F) != 1:
            raise Plasmid_Exception('There should be exactly one forward {0} site in {1} but {2} were found!' % (restriction_enzyme, part_plasmid.plasmid_name or part_plasmid.get_id(), sequence.count(site_F)))
        if sequence.count(site_R) != 1:
            raise Plasmid_Exception('There should be exactly one reverse {0} site in {1} but {2} were found!' % (restriction_enzyme, part_plasmid.plasmid_name or part_plasmid.get_id(), sequence.count(site_R)))

        left_overhangs = set([p['Part_Type_overhang_5'] for p in plasmid_parts])
        right_overhangs = set([p['Part_Type_overhang_3'] for p in plasmid_parts])

        assert(len(left_overhangs) == len(right_overhangs))
        leftoverhangs = left_overhangs.symmetric_difference(right_overhangs) # Leftover overhangs? Get it?

        left_overhang, right_overhang = None, None
        for p in plasmid_parts:
            if p['Part_Type_overhang_5'] in leftoverhangs:
                left_overhang = p['Part_Type_overhang_5']
            if p['Part_Type_overhang_3'] in leftoverhangs:
                right_overhang = p['Part_Type_overhang_3']

        site_F_position = sequence.find(left_overhang, sequence.find(site_F))
        site_R_position = sequence.find(right_overhang, sequence.find(site_R) - 10) + 4

        # Circular permutation so that forward cut site is always upstream of reverse cut site in linear sequence
        if site_F_position > site_R_position:
            sequence = sequence[site_R_position + 8:] + sequence[:site_R_position + 8]

        new_part_sequence = sequence[sequence.find(left_overhang, sequence.find(site_F)):sequence.find(right_overhang,sequence.find(site_R) - 10) + 4]
        if new_part_sequence not in part_sequences:
            part_sequences.append(new_part_sequence)

        return part_sequences, part_types



    def fetch_cassette_parts(self, tsession, input_sequences, table_info=None):
        '''
        Get sequences and part numbers for listed part plasmids

        todo: This is an old version of the function that I left here so as not to break the cassette assembly code.
              It should be either merged with fetch_plasmid_parts or changed to take some of the speed-ups in that function
              or maybe this code belongs in Plasmid_Utilities (there is already a version of it there so code drift is a
              problem).

        :param tsession:
        :param input_sequences:
        :param table_info:
        :return:
        '''

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
        part_types = []

        for plasmid, part_plasmid, part_plasmid_part, part_type in part_plasmids_query:

            #print "%s\t%s\t%s\t%s" % (plasmid.plasmid_name, plasmid.creator, plasmid.creator_entry_number, part_plasmid_part.part_number)

            part_types.append(part_type.part_number)

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

            new_part_sequence = sequence_upper[sequence_upper.find(left_overhang, sequence_upper.find(site_F)):sequence_upper.find(right_overhang,sequence_upper.find(site_R) - 10) + 4]
            if new_part_sequence not in part_sequences:
                part_sequences.append(new_part_sequence)

        if table_info:
            return table_info
        else:
            return part_sequences, part_types # todo: this is different from James's function


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


    def get_details(self, tsession, engine, cache = None, fresh = False):

        if fresh or '__details__' not in self.__dict__:

            color = tsession.query(Feature_Type).filter(Feature_Type.Feature_type == self.Feature_type).one().color

            protein_sequence_translation = None
            try:
                if len(self.Feature_sequence) % 3 == 0:
                    protein_sequence_translation = translate_codons(self.Feature_sequence)
            except: pass

            self.__dict__['__details__'] = dict(
                    name = self.Feature_name,
                    hash = self.MD5_hash,
                    type = self.Feature_type,
                    sequence = self.Feature_sequence,
                    protein_sequence = protein_sequence_translation,
                    description = self.description,
                    color = color,
            )

        return self.__dict__['__details__']


class Feature_Type(DeclarativeBasePlasmid):
    __tablename__ = 'Feature_Type'

    Feature_type = Column(Unicode(200, collation="utf8_bin"), nullable=False, primary_key=True)
    color = Column(Unicode(7, collation="utf8_bin"), nullable=False)
    designable = Column(Integer, nullable=False)


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


class Other_Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Other_Plasmid'

    _required_fields = ['creator', 'creator_entry_number', 'resistance']

    creator_entry_number = Column(Integer, ForeignKey('Plasmid.creator_entry_number'), nullable=False, primary_key=True)
    creator = Column(Unicode(5), ForeignKey('Plasmid.creator'), nullable=False, primary_key=True)
    resistance = Column(Unicode(100, collation="utf8_bin"), nullable=False)
    vector = Column(Unicode(32, collation = "utf8_bin"), nullable = True)

    plasmid = relationship('Plasmid', primaryjoin='and_(Plasmid.creator == Other_Plasmid.creator, Plasmid.creator_entry_number == Other_Plasmid.creator_entry_number)')


    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Other_Plasmid(**input_dict)

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
        d['__id__'] = self.plasmid.get_id() # this is easier to work with for the website hashtables
        return d


    def __repr__(self):
        return 'Internal numbering: {0} {1}'.format(self.creator, self.creator_entry_number)


class Design_Plasmid(DeclarativeBasePlasmid):
    __tablename__ = 'Design_Plasmid'

    _required_fields = ['creator', 'creator_entry_number', 'resistance']

    creator_entry_number = Column(Integer, ForeignKey('Plasmid.creator_entry_number'), nullable = False, primary_key = True)
    creator = Column(Unicode(5), ForeignKey('Plasmid.creator'), nullable = False, primary_key = True)
    resistance = Column(Unicode(100, collation = "utf8_bin"), nullable = False)
    vector = Column(Unicode(32, collation = "utf8_bin"), nullable = True)
    is_template = Column(TINYINT(1), default = 0, nullable = False)

    plasmid = relationship('Plasmid', primaryjoin = 'and_(Plasmid.creator == Design_Plasmid.creator, Plasmid.creator_entry_number == Design_Plasmid.creator_entry_number)')


    @staticmethod
    def add(tsession, input_dict, silent = True):
        try:
            db_record_object = Design_Plasmid(**input_dict)

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
        d['__id__'] = self.plasmid.get_id()  # this is easier to work with for the website hashtables
        return d


    def __repr__(self):
        return 'Internal numbering: {0} {1}'.format(self.creator, self.creator_entry_number)


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


    def get_details(self, tsession, engine):
        '''Returns details about a feature.'''
        feature_details = self.feature.get_details(tsession, engine)
        plasmid_sequence = self.plasmid.sequence.upper().strip()
        feature_sequence = self.feature.Feature_sequence.upper().strip()
        feature_type = self.feature.ftype
        indices = []
        feature_details['indices'] = indices

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


    def get_details(self, tsession, engine):
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
    mutation = Column(String(3), nullable=False)
    position = Column(Integer, nullable=False)
    wt_AA = Column(String(1), nullable=False)
    mut_AA = Column(String(1), nullable=False)
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


class Publication(DeclarativeBasePlasmid):
    __tablename__ = 'Publication'

    ID = Column(Integer, nullable = False, primary_key = True)
    Title = Column(Unicode(256), nullable = True)
    Publication = Column(String(256), nullable = True)
    Volume = Column(String(8), nullable = True)
    Issue = Column(String(8), nullable = True)
    StartPage = Column(String(16), nullable = True)
    EndPage = Column(String(16), nullable = True)
    PublicationYear = Column(Integer, nullable = True)
    PublicationDate = Column(DateTime, nullable = True)
    DOI = Column(String(64), nullable=True)
    URL = Column(String(128), nullable=True)
    ISSN = Column(Integer, nullable=True)
    PubMedID = Column(String(64), nullable=True)
    RecordType = Column(String(128), nullable=True)
    Notes = Column(Text, nullable=True)
    RIS = Column(Text, nullable=True)

    authors = relationship('PublicationAuthor', primaryjoin='PublicationAuthor.PublicationID == Publication.ID', order_by=lambda: PublicationAuthor.AuthorOrder)


    def get_details(self):
        '''This function is used by the web interface.'''
        d = row_to_dict(self)
        d['authors'] = [dict(
                FirstName = a.FirstName,
                MiddleNames = a.MiddleNames,
                Surname = a.Surname,
            ) for a in self.authors]
        return d


    @staticmethod
    def add(tsession, publication_record, creator, creator_entry_number, description, silent=True, allow_existing_publication_record=False):

        try:
            # Add the Publication record
            d = {}
            expected_fields = {}
            optional_fields = {
                'Title' : 'Title',
                'PublicationName' : 'Publication',
                'Volume' : 'Volume',
                'Issue' : 'Issue',
                'StartPage': 'StartPage',
                'EndPage'  : 'EndPage',
                'PublicationYear' : 'PublicationYear',
                'PublicationDate' : 'PublicationDate',
                'DOI': 'DOI',
                'URL': 'URL',
                'ISSN': 'ISSN',
                'PubMedID': 'PubMedID',
                'RecordType': 'RecordType',
                'Notes': 'Notes',
                'RIS': 'RIS',
            }
            for k in expected_fields:
                assert(k in publication_record and publication_record[k] != None)
                d[expected_fields[k]] = publication_record[k]
            for k in optional_fields:
                d[optional_fields[k]] = publication_record.get(k)
            new_publication = Publication(**d)

            existing_hits = []
            if new_publication.DOI != None:
                existing_hits.extend([r for r in tsession.query(Publication).filter(Publication.DOI == new_publication.DOI)])
            if new_publication.URL != None:
                existing_hits.extend([r for r in tsession.query(Publication).filter(Publication.URL == new_publication.URL)])
            if new_publication.ISSN != None:
                existing_hits.extend([r for r in tsession.query(Publication).filter(Publication.ISSN == new_publication.ISSN)])
            if new_publication.PubMedID != None:
                existing_hits.extend([r for r in tsession.query(Publication).filter(Publication.PubMedID == new_publication.PubMedID)])
            print(existing_hits)
            if existing_hits:
                if not allow_existing_publication_record:
                    raise DuplicatePublicationException()
                else:
                    new_publication = existing_hits[0]
            else:
                tsession.add(new_publication)
                tsession.flush()

                # Add the PublicationAuthor records
                authors = publication_record.get('authors', [])
                assert(len(authors) > 0)
                for a in authors:
                    a['PublicationID'] = new_publication.ID
                    new_publication_author = PublicationAuthor(**a)
                    tsession.add(new_publication_author)
                    tsession.flush()

            # Add the Publication_Plasmid record
            for r in tsession.query(Publication_Plasmid).filter(and_(
                    Publication_Plasmid.PublicationID == new_publication.ID,
                    Publication_Plasmid.creator == creator,
                    Publication_Plasmid.creator_entry_number == creator_entry_number,
                )):
                r.Deprecated = 1

            new_plasmid_publication = Publication_Plasmid(**dict(
                PublicationID = new_publication.ID,
                creator = creator,
                creator_entry_number = creator_entry_number,
                Description = description,
                Deprecated = 0,
            ))
            tsession.add(new_plasmid_publication)
            tsession.flush()

            return new_publication

        except:
            raise


class PublicationAuthor(DeclarativeBasePlasmid):
    __tablename__ = 'PublicationAuthor'

    PublicationID = Column(Integer, ForeignKey('Publication.ID'), nullable = False, primary_key = True)
    AuthorOrder = Column(Integer, nullable = False, primary_key = True)
    FirstName = Column(Unicode(64), nullable = False)
    MiddleNames = Column(Unicode(64), nullable = True)
    Surname = Column(Unicode(64), nullable = True)


class Plasmid_Feature_Design(DeclarativeBasePlasmid):
    __tablename__ = 'Plasmid_Feature_Design'

    parent_creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Users.ID'), nullable=False, primary_key=True)
    parent_creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    child_creator = Column(Unicode(5, collation="utf8_bin"), ForeignKey('Users.ID'), nullable=False, primary_key=True)
    child_creator_entry_number = Column(Integer, nullable=False, primary_key=True)
    feature_ID = Column(Integer, nullable=False, primary_key=True)

    @staticmethod
    def add(tsession, input_dict, silent=True):
        try:
            db_record_object = Plasmid_Feature_Design(**input_dict)

            if not silent:
                colortext.pcyan('Adding this record:')
                print(db_record_object)
                print('')

            tsession.add(db_record_object)
            tsession.flush()

            return db_record_object
        except:
            raise


class Plasmid_Type(DeclarativeBasePlasmid):
    __tablename__ = 'Plasmid_Type'

    plasmid_type = Column(Unicode(64), nullable=False, primary_key=True)
    show_on_add_plasmid_page = Column(TINYINT(1), default = 0, nullable=False)