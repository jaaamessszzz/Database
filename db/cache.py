import traceback

from sqlalchemy.sql import func, desc, join
from sqlalchemy import and_, or_, func

from model import *

class CacheMissException(Exception): pass

class Cache(object):


    def __init__(self, tsession, engine):

        self.tsession = tsession
        self.engine = engine

        self.features = {}
        self.plasmids = {}
        self.plasmid_order = []


    def load_features(self, tsession = None, engine = None):

        engine = engine or self.engine
        tsession = tsession or self.tsession

        for f in tsession.query(Feature):
            self.features[f.Feature_name] = f


    def load_plasmid(creator, creator_entry_number, tsession = None, engine = None):
        engine = engine or self.engine
        tsession = tsession or self.tsession

        if k not in self.plasmids:
            k = (creator, creator_entry_number)
            p = tsession.query(Plasmid).filter(Plasmid.creator == creator, Plasmid.creator_entry_number == creator_entry_number).one()
            self.plasmid_order.append(k)
            self.plasmids[k] = p
        return self.plasmids[k]


    def load_plasmids(self, tsession = None, engine = None):

        engine = engine or self.engine
        tsession = tsession or self.tsession

        for p in tsession.query(Plasmid).order_by(desc(Plasmid.date), desc(Plasmid.creator_entry_number)): # todo: parameterize sorting
            k = (p.creator, p.creator_entry_number)
            self.plasmid_order.append(k)
            self.plasmids[k] = p


    def counts(self):
        return dict(
            features = len(self.features),
            plasmids = len(self.plasmids),
        )


    def get_feature(self, feature_name):
        try: return self.features[feature_name]
        except: raise CacheMissException()


    def get_plasmid(self, creator, creator_entry_number, tsession = None, engine = None):
        try:
            return self.plasmids[(creator, creator_entry_number)]
        except Exception, e:
            pass
            #print(str(e))
            #print(traceback.format_exc())

            try:
                return self.load_plasmid(creator, creator_entry_number, tsession = tsession, engine = engine)
            except Exception, e:
                print(str(e))
                print(traceback.format_exc())
                raise CacheMissException()


    def get_plasmid_tpl(self, creator, creator_entry_number):
        try: return self.plasmids[(creator, creator_entry_number)]
        except: raise CacheMissException()

