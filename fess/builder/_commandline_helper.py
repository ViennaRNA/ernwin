from __future__ import absolute_import
import logging
log=logging.getLogger(__name__)


def replica_substring(string, replica_nr):
    log.debug("String is %r", string)
    if replica_nr is None or '@' not in string:
        return string
    else:
        try:
            return string.split('@')[replica_nr]
        except IndexError:
            raise IndexError("Not enough '@'-characters in string '{}'."
                             "For replica-exchange, the number of '@'-separated "
                             "contributions has to be equal to the number "
                             "of replicas.".format(string))
