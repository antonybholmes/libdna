import libdna
import logging
import unittest
import sys

class TestDecode(unittest.TestCase):
    def decode(self):
        logging.basicConfig(stream=sys.stderr)
        log = logging.getLogger(__name__).setLevel( logging.DEBUG )

        location = "chr3:187446721-187447977"
    
        dna = libdna.DNA2Bit('/ifs/scratch/cancer/Lab_RDF/abh2138/references/ucsc/assembly/2bit')
    
        #dna.fasta(location)
    
        dna.fasta(('chr1', 14453, 14453 + 350))
        
        s = dna.dna('chr1:14453-15000')
        
        log.debug(s)
        log.debug('what')
        
        self.assertTrue(isinstance(s, str))
