import os
from abc import ABC, abstractmethod
import libdna

class DNA(ABC):
    @abstractmethod
    def dna(self, *args):
        raise NotImplementedError

    def fasta(self, loc, mask='mask'):
        """
        Prints a fasta representation of a sequence.
        
        Parameters
        ----------
        l : tuple (str, int, int)
            location chr, start, and end
        mask : str, optional
            Either 'upper', 'lower', or 'n'. If 'lower', poor quality bases 
            will be converted to lowercase.
        """
        
        l = libdna.parse_loc(loc)
        
        print('>{}'.format(l))
        print(self.dna(l, mask=mask))
        
        
class DNAStr(DNA):
    def __init__(self, dir):
        self.__dir = dir
    
    def dna(self, loc):
        l = libdna.parse_loc(loc)
        
        sstart = l.start - 1
        send = l.end - 1
      
        file = os.path.join(self.dir, l.chr + ".txt")
        
        f = open(file, 'r')
      
        f.seek(sstart)
      
        length = send - sstart + 1
        
        seq = f.read(length)
    
        f.close()
      
        return seq
    
    
class DNA2Bit(DNA):
    def __init__(self, dir):
        self.__dir = dir
        self.__data = []
        self.__file = ''
        self.__n_data = []
        self.__n_file = ''
        self.__mask_data = []
        self.__mask_file = ''
        
    @staticmethod
    def decoded_base(base):
        if base == 0:
            return 'A'
        elif base == 1:
            return 'C'
        elif base == 2:
            return 'G'
        else:
            return 'T'
        
    @staticmethod
    def _read1bit(d, l):
        """
        Read data from a 1 bit file where each byte encodes 8 bases.
        
        Parameters
        ----------
        file : str
            A binary file to read.
        l : tuple
            chr, start, end
        
        Returns
        -------
        list
            list of 1s and 0s of length equal to the number of bases in
            the location.
        """

        s = l.start - 1
        si = s // 8

        length = l.end - l.start + 1

        ret = []
        
        b = s
        bi = si
        
        for i in range(0, length):
            block = b % 8
            
            if block == 0:
                v = (d[bi] >> 7)
            elif block == 1:
                v = (d[bi] >> 6)
            elif block == 2:
                v = (d[bi] >> 5)
            elif block == 3:
                v = (d[bi] >> 4)
            elif block == 4:
                v = (d[bi] >> 3)
            elif block == 5:
                v = (d[bi] >> 2)
            elif block == 6:
                v = (d[bi] >> 1)
            else:
                v = d[bi]
                bi += 1
            
            # Only care about the lowest bit
            v &= 1
            
            ret.append(v)
            
            b += 1
    
        return ret
    
    def _read_dna(self, l):
        """
        Read DNA from a 2bit file where each base is encoded in 2bit 
        (4 bases per byte).
        
        Parameters
        ----------
        l : tuple
            Location tuple
        
        Returns
        -------
        list
            Array of base chars
        """
        
        file = os.path.join(self.__dir, l.chr + ".dna.2bit")
        
        if not os.path.exists(file):
            return []

        if file != self.__file:
            print('Caching {}...'.format(file))
            f = open(file, 'rb')
            self.__data = f.read()
            f.close()
            self.__file = file

        s = l.start - 1
        si = s // 4
        #ei = (end - 1) // 4

        length = l.end - l.start + 1
    
        ret = []
        
        b = s
        bi = si
        
        for i in range(0, length):
            block = b % 4
            
            if block == 0:
                v = (self.__data[bi] >> 6)
            elif block == 1:
                v = (self.__data[bi] >> 4)
            elif block == 2:
                v = (self.__data[bi] >> 2)
            else:
                v = self.__data[bi]
                
                # Reached end of byte so we are moving into the next byte
                bi += 1
            
            # Only care about the lowest 2 bits
            v &= 3
            
            base = DNA2Bit.decoded_base(v)
        
            ret.append(base)
            
            b += 1
    
        return ret
    
    def _read_n(self, l, ret):
        """
        Reads 'N' mask from 1 bit file to convert bases to 'N'. In the
        2 bit file, 'N' or any other invalid base is written as 'A'.
        Therefore the 'N' mask file is required to correctly identify where
        invalid bases are.
        
        Parameters
        ----------
        l : tuple
            location
        ret : list
            List of bases which will be modified in place.
        """
        
        file = os.path.join(self.__dir, l.chr + ".n.1bit")
        
        if not os.path.exists(file):
            return
        
        if file != self.__n_file:
            print('Caching {}...'.format(file))
            f = open(file, 'rb')
            self.__n_data = f.read()
            f.close()
            self.__n_file = file
        
        d = DNA2Bit._read1bit(self.__n_data, l)
        
        for i in range(0, len(ret)):
            if d[i] == 1:
                ret[i] = 'N'
                    
    def _read_mask(self, l, ret, mask='upper'):
        """
        Reads mask from 1 bit file to convert bases to identify poor quality
        bases that will either be converted to lowercase or 'N'. In the
        2 bit file, 'N' or any other invalid base is written as 'A'.
        Therefore the 'N' mask file is required to correctly identify where
        invalid bases are.
        
        Parameters
        ----------
        l : tuple
            location
        ret : list
            list of bases which will be modified in place
        mask : str, optional
            Either 'upper', 'lower', or 'n'. If 'lower', poor quality bases 
            will be converted to lowercase.
        """
        
        if mask == 'upper':
            return
         
        file = os.path.join(self.__dir, l.chr + ".mask.1bit")
             
        if not os.path.exists(file):
            return
        
        if file != self.__mask_file:
            print('Caching {}...'.format(file))
            f = open(file, 'rb')
            self.__mask_data = f.read()
            f.close()
            self.__mask_file = file
        
        d = DNA2Bit._read1bit(self.__mask_data, l)
        
        if mask == 'lower':
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = ret[i].lower()
        else:
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = 'N'
    
    def dna(self, loc, mask='upper'):
        """
        Returns the DNA for a location.
        
        Parameters
        ----------
        mask : str, optional
            Indicate whether masked bases should be represented as is
            ('upper'), lowercase ('lower'), or as N ('n')
        
        Returns
        -------
        list
            List of base chars.
        """
        
        l = libdna.parse_loc(loc)
            
        ret = self._read_dna(l)
        
        self._read_n(l, ret)
            
        self._read_mask(l, ret, mask=mask)
      
        return ''.join(ret)

    
    def merge_read_pair_seq(self, r1, r2):
        """
        Merge the sequence of two reads into one continuous read either
        by inserting the missing DNA, or joining on the common sequence.
        
        Parameters
        ----------
        r1 : libsam.Read
            Read 1
        r2 : libsam.Read
            Read 2
        """
        
        s1 = r1.pos # + 1
            
        # end of first read
        e1 = s1 + r1.length - 1
    
        # start of second read
        s2 = r2.pos # + 1
        
        e2 = s2 + r2.length - 1
    
        inner = s2 - e1 - 1
    
        if inner >= 0:
            seq = self.dna((r1.chr, s1, e2))
        else:
            # Reads overlap so concatenate the first read with the
            # portion of the second read that is not overlapping
            # (inner is negative so flip sign for array indexing)
            seq = r1.seq + r2.seq[-inner:]
            
        return seq
