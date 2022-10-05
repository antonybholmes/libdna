import os
from abc import ABC, abstractmethod
import libdna
import sys
import s3fs

# Use ord('A') etc to get ascii values
DNA_UC_DECODE_DICT = {0:65, 1:67, 2:71, 3:84}
DNA_LC_DECODE_DICT = {0:97, 1:99, 2:103, 3:116}
DNA_UC_TO_LC_MAP = {65:97, 67:99, 71:103, 84:116, 78:110}
DNA_N_UC = 78
DNA_N_LC = 110
DNA_COMP_DICT = {65:84, 67:71, 84:65, 71:67, 97:116, 99:103, 116:97, 103:99, 78:78, 110:110}
EMPTY_BYTEARRAY = bytearray(0)

class DNA(ABC):
    @abstractmethod
    def dna(self, *args):
        raise NotImplementedError

    def fasta(self, loc, mask='upper'):
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
        
    def to_fasta(self, loc, file, mask='upper'):
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
        
        f = open(file, 'w')
        f.write('>{}\n'.format(l))
        f.write('{}\n'.format(self.dna(l, mask=mask)))
        f.close()
        
        
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
        
    @property
    def dir(self):
        return self.__dir
    
    
    @staticmethod
    def rev_comp(dna):
        """
        Parameters
        ----------
        dna : bytearray
            dna sequence to be reverse complemented
        """
        
        i2 = len(dna) - 1
        
        l = len(dna) // 2
        
        for i in range(0, l):
            b = DNA_COMP_DICT[dna[i]]
            dna[i] = DNA_COMP_DICT[dna[i2]]
            dna[i2] = b
            i2 -= 1
            
            
    def _read1bit(self, d, loc, offset=False):
        """
        Read data from a 1 bit file where each byte encodes 8 bases.
        
        Parameters
        ----------
        d : array
            byte array
        l : tuple
            chr, start, end
        
        Returns
        -------
        list
            list of 1s and 0s of length equal to the number of bases in
            the location.
        """
        
        if d is None:
            return []

        s = loc.start - 1

        length = loc.end - loc.start + 1

        ret = [0] * length
        
        if offset:
            bi = s // 8
        else:
            bi = 0
            
        for i in range(0, length):
            block = s % 8
            
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
            
            ret[i] = v
            
            s += 1
    
        return ret
    
    
    def _read2bit(self, d, loc, offset=False):
        """
        Read DNA from a 2bit file where each base is encoded in 2bit 
        (4 bases per byte).
        
        Parameters
        ----------
        d: bytes array
             Encoded dna.
        loc : tuple
            Location (chr, start, end)
        
        Returns
        -------
        list
            Array of base chars
        """
        
        if d is None:
            return EMPTY_BYTEARRAY
        
        #print(d, loc)
        
        s = loc.start - 1
        
        ret = bytearray([0] * loc.length) #[]
        
        if offset:
            bi = s // 4
        else:
            bi = 0
        
        for i in range(0, loc.length):
            block = s % 4
            
            if block == 0:
                v = (d[bi] >> 6)
            elif block == 1:
                v = (d[bi] >> 4)
            elif block == 2:
                v = (d[bi] >> 2)
            else:
                v = d[bi]
                
                # Reached end of byte so we are moving into the next byte
                bi += 1
            
            # Only care about the lowest 2 bits
            v &= 3
            
            ret[i] = DNA_UC_DECODE_DICT[v]
                
            s += 1

  
        return ret
    
    
    def _read_dna(self, loc, lowercase=False):
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
        
        file = '{}.dna.2bit'.format(loc.chr)
        
        s = loc.start - 1
        e = s + loc.length
        bs = s // 4
        be = e // 4
        l = be - bs + 1
        
#        f = open(file, 'rb')
#        f.seek(bs) #(loc.start - 1) // 4)
#        # read bytes into buffer
#        data = f.read(l) #loc.length // 4 + 2)
#        f.close()
        
        data = self.read_data(file, bs, l)
        
        print(loc, data, file=sys.stderr)
        
        return self._read2bit(data, loc)
    
    
    def _read_1bit_file(self, file, loc):
        """
        Load data from 1 bit file into array
        
        Parameters
        ----------
        file : str
            1bit filename
        loc : libdna.Loc
            dna location
        
        Returns
        -------
        bytes
            byte array from file where each byte represents 8 bases.
        """
        
        s = loc.start - 1
        e = s + loc.length
        bs = s // 8
        be = e // 8
        n = be - bs + 1
        
        
        data = self.read_data(file, bs, n)
        
#        f = open(file, 'rb')
#        f.seek(bs)  #(l.start - 1) // 8)
#        # read length + 2 because we need the extra byte in case the start
#        # position lies mid way through a byte. Imagine a sequence 10 bp long
#        # starting at position 4. 10 // 8 + 1 = 2 bytes of data required to
#        # store this. Since we pick the closest byte as the start, this will
#        # be position 0 (3 // 8 = 0). The length is 2 bytes spanning bytes 0
#        # and 1, but because of the start at 3, our sequence spans byte 3 so
#        # we need to buffer an extra byte for cases where the start does not
#        # match the start of a byte
#        data = f.read(l) #l.length // 8 + 2)
#        f.close()
        return data
    
    
    def _read_n(self, loc, ret):
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
        
        file = '{}.n.1bit'.format(loc.chr)
        
        data = self._read_1bit_file(file, loc)
        
        d = self._read1bit(data, loc)
        
        for i in range(0, len(ret)):
            if d[i] == 1:
                ret[i] = DNA_N_UC #'N'
                
                
    def read_data(self, file, seek, n):
        """
        Reads data from a file source
        
        Parameter
        ---------
        file : str
            Relative path to file
        seek : int
            Start offset in bytes
        n : int
            Amoount of data to read in bytes
        
        Returns
        -------
        bytearray
            Data from file
        """
        file = os.path.join(self.dir, file).lower()
        
        if not os.path.exists(file):
            return None
        
        f = open(file, 'rb')
        f.seek(seek)
        data = f.read(n) #l.length // 8 + 2)
        f.close()
        return data
    
        
    def _read_mask(self, loc, ret, mask='upper'):
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
        
        if mask.startswith('u'):
            return
         
        file = '{}.mask.1bit'.format(loc.chr)
        
        data = self._read_1bit_file(file, loc)
        
        d = self._read1bit(data, loc)
        
        if mask.startswith('l'):
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = DNA_UC_TO_LC_MAP[ret[i]] #ret[i].lower()
        else:
            # Use N as mask
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = DNA_N_UC #'N'
                    
    
    def dna(self, loc, mask='lower', rev_comp=False, lowercase=False):
        """
        Returns the DNA for a location.
        
        Parameters
        ----------
        loc : str or tuple
            Genomic Location
        mask : str, optional
            Indicate whether masked bases should be represented as is
            ('upper'), lowercase ('lower'), or as N ('n')
        lowercase : bool, optional
            Indicates whether sequence should be displayed as upper or
            lowercase. Default is False so sequence is uppercase. Note that
            this only affects the reference DNA and does not affect the
            mask.
        
        Returns
        -------
        list
            List of base chars.
        """
        
        l = libdna.parse_loc(loc)
            
        ret = self._read_dna(l, lowercase=lowercase)
        
        self._read_n(l, ret)
        
        self._read_mask(l, ret, mask=mask)
        
        if rev_comp:
            DNA2Bit._rev_comp(ret)

        ret = ret.decode('utf-8')
        
        if lowercase:
            ret = ret.lower()
      
        return ret

    
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
            seq = self.dna(libdna.Loc(r1.chr, s1, e2))
        else:
            # Reads overlap so concatenate the first read with the
            # portion of the second read that is not overlapping
            # (inner is negative so flip sign for array indexing)
            seq = r1.seq + r2.seq[-inner:]
            
        return seq

  
class S3DNA2Bit(DNA2Bit):
    def __init__(self, bucket, dir):
        super().__init__(dir)
        self.__bucket = bucket
        self.__fs = s3fs.S3FileSystem(anon=True)
        
    @property
    def bucket(self):
        return self.__bucket
    
    def read_data(self, file, seek, n):
        file = os.path.join(self.__bucket, self.dir, file).lower()
        
        f = self.__fs.open(file, 'rb')
        f.seek(seek)
        data = f.read(n)
        return data
    
    
class CachedDNA2Bit(DNA2Bit):
    def __init__(self, dir):
        super().__init__(dir)
        
        self.__data = []
        self.__file = ''
        self.__n_data = []
        self.__n_file = ''
        self.__mask_data = []
        self.__mask_file = ''
        
    
    def _read_dna(self, loc, lowercase=False):
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
        
        file = os.path.join(self.dir, loc.chr + ".dna.2bit")
        
        if not os.path.exists(file):
            print(file, 'does not exist.')
            return EMPTY_BYTEARRAY

        if file != self.__file:
            print('Caching {}...'.format(file))
            self.__file = file
            # Load file into memory
            f = open(file, 'rb')
            self.__data = f.read()
            f.close()
            
            
        return self._read2bit(self.__data, loc, offset=True)
    
    
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
        
        file = os.path.join(self.dir, l.chr + ".n.1bit")
        
        if not os.path.exists(file):
            return
        
        if file != self.__n_file:
            print('Caching {}...'.format(file))
            f = open(file, 'rb')
            self.__n_data = f.read()
            f.close()
            self.__n_file = file
        
        d = self._read1bit(self.__n_data, l, offset=True)
        
        for i in range(0, len(ret)):
            if d[i] == 1:
                ret[i] = DNA_N_UC #'N'
                
                    
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
        
        if mask.startswith('u'):
            return
         
        file = os.path.join(self.dir, l.chr + ".mask.1bit")
             
        if not os.path.exists(file):
            return
        
        if file != self.__mask_file:
            print('Caching {}...'.format(file))
            f = open(file, 'rb')
            self.__mask_data = f.read()
            f.close()
            self.__mask_file = file
        
        d = self._read1bit(self.__mask_data, l, offset=True)
        
        if mask.startswith('l'):
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = DNA_UC_TO_LC_MAP[ret[i]] #ret[i].lower()
        else:
            # Use N as mask
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = DNA_N_UC #'N'



