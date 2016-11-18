if __name__ == '__main__':
    import sys
    sys.path.insert(0, '..')

import traceback

from Utilities.Plasmid_Utilities import Plasmid_Utilities
from db.model import Plasmid

def cincinnati_squish_em_all():
    '''Note: nuke_plasmids is a very dangerous function to call.'''
    try:
        utilities = Plasmid_Utilities()
        baes = [(bae.creator, bae.creator_entry_number) for bae in utilities.tsession.query(Plasmid).filter(Plasmid.creator == u'PARP')]
        for bae in baes:
            assert(bae[0] == 'PARP')
        baes = [('PARP', 77)]
        utilities.nuke_plasmids(baes)
        utilities.tsession.close()
    except Exception, e:
        if utilities and utilities.tsession:
            utilities.tsession.close()
        print(str(e))
        print(traceback.format_exc())


if __name__ == '__main__':
    cincinnati_squish_em_all()

