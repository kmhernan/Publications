#!/usr/bin/env python
# Although not required in any sense, share the love and pass on attribution
# when using or modifying this code.
#
# To the extent possible under law, the author(s) have dedicated all copyright
# and related and neighboring rights to this software to the public domain
# worldwide. This software is distributed without any warranty.
#
# You should have received a copy of the CC0 Public Domain Dedication along with
# this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>

import gzip

class Reference(object):
    """
    This class represents a FASTA formatted reference file.
    The constructor takes a filehandle and creates a generator
    which returns a subsequent scaffold.
    """
    def __init__(self, handle, compressed = False):
        self.handle = handle
        self.compressed = compressed

    def __iter__(self):
        scaff = ''
        seq   = []
        if self.compressed:
            for line in gzip.open(self.handle, 'rb'):
                if line.startswith('>') and not scaff:
                    scaff = line.rstrip().replace('>', '')
                elif line.startswith('>') and scaff:
                    self.scaffold = scaff
                    self.sequence = ''.join(seq)
                    yield self

                    scaff = line.rstrip().replace('>', '')
                    seq   = []
                else:
                    seq.append(line.rstrip())
            self.scaffold = scaff
            self.sequence = ''.join(seq)
            yield self

        else:
            for line in open(self.handle, 'rU'):
                if line.startswith('>') and not scaff:
                    scaff = line.rstrip().replace('>', '')
                elif line.startswith('>') and scaff:
                    self.scaffold = scaff
                    self.sequence = ''.join(seq)
                    yield self

                    scaff = line.rstrip().replace('>', '')
                    seq   = []
                else:
                    seq.append(line.rstrip())
            self.scaffold = scaff
            self.sequence = ''.join(seq)
            yield self

    def len(self):
        return len(self.sequence)

    def write(self, ostream):
        ostream.write('>{0}\n{1}\n'.format(self.scaffold, self.sequence))
