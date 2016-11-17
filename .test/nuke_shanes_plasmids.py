if __name__ == '__main__':
    import sys
    sys.path.insert(0, '..')

from Utilities.Plasmid_Utilities import Plasmid_Utilities
from Utilities.Plasmid_View_Tools import Plasmid_View_Tools
from db.interface import DatabaseInterface
from db.model import Feature, Plasmid_File, Plasmid, Part_Plasmid_Part, CDS_Mutant_Constituent
from sqlalchemy import and_

import sys
import pprint
import traceback

def cincinnati_squish_em_all():
    '''Note: nuke_plasmids is a very dangerous function to call.'''

    # Create up the database session
    try:
        utilities = Plasmid_Utilities()
        baes = [(bae.creator, bae.creator_entry_number) for bae in utilities.tsession.query(Plasmid).filter(Plasmid.creator == u'PARP')]
        for bae in baes:
            assert(bae[0] == 'PARP')
        #utilities.nuke_plasmids(baes)
        utilities.tsession.close()
    except Exception, e:
        if utilities and utilities.tsession:
            utilities.tsession.close()
        print(str(e))
        print(traceback.format_exc())


if __name__ == '__main__':
    import sys
    sys.path.insert(0, '..')

    cincinnati_squish_em_all()

