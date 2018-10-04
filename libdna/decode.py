import os
from abc import ABC, abstractmethod
import libdna

# Use ord('A') etc to get ascii values
DNA_UC_DECODE_DICT = {0:65, 1:67, 2:71, 3:84}
DNA_LC_DECODE_DICT = {0:97, 1:99, 2:103, 3:116}
DNA_UC_TO_LC_MAP = {65:97, 67:99, 71:103, 84:116, 78:110}
DNA_N_UC = 78
DNA_N_LC = 110
DNA_COMP_DICT = {65:84, 67:71, 84:65, 71:67, 97:116, 99:103, 116:97, 103:99, 78:78, 110:110}

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
            
        
    @staticmethod
    def _read1bit(d, l, offset=False):
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

        s = l.start - 1

        length = l.end - l.start + 1

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
    
    
    @staticmethod
    def _read2bit(d, l, offset=False):
        """
        Read DNA from a 2bit file where each base is encoded in 2bit 
        (4 bases per byte).
        
        Parameters
        ----------
        d:
        l : tuple
            Location tuple
        
        Returns
        -------
        list
            Array of base chars
        """
        
        s = l.start - 1
        
        ret = bytearray([0] * l.length) #[]
        
        if offset:
            bi = s // 4
        else:
            bi = 0
        
        for i in range(0, l.length):
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
    
        print('dfdf', type(ret))
        return ret
    
    
    def _read_dna(self, l, lowercase=False):
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
        
        file = os.path.join(self.dir, l.chr + ".dna.2bit")
        
        print(file)
        
        if not os.path.exists(file):
            return bytearray([])
       
        f = open(file, 'rb')
        f.seek((l.start - 1) // 4)
        # read bytes into buffer
        data = f.read(l.length // 4 + 1)
        f.close()
        
        return DNA2Bit._read2bit(data, l)
    
    
    @staticmethod
    def _read_1bit_file(file, l):
        """
        Load data from 1 bit file into array
        
        Parameters
        ----------
        file : str
            1bit filename
        l : libdna.Loc
            dna location
        
        Returns
        -------
        bytes
            byte array from file where each byte represents 8 bases.
        """
        
        f = open(file, 'rb')
        f.seek((l.start - 1) // 8)
        data = f.read(l.length // 8 + 1)
        f.close()
        return data
    
    
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
        
        data = DNA2Bit._read_1bit_file(file, l)
        
        d = DNA2Bit._read1bit(data, l)
        
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
         
        file = os.path.join(self.__dir, l.chr + ".mask.1bit")
             
        if not os.path.exists(file):
            return
        
        data = DNA2Bit._read_1bit_file(file, l)
        
        d = DNA2Bit._read1bit(data, l)
        
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
        
        print('f', type(ret))
        
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
            seq = self.dna((r1.chr, s1, e2))
        else:
            # Reads overlap so concatenate the first read with the
            # portion of the second read that is not overlapping
            # (inner is negative so flip sign for array indexing)
            seq = r1.seq + r2.seq[-inner:]
            
        return seq
    
    
class CachedDNA2Bit(DNA2Bit):
    def __init__(self, dir):
        super().__init__(dir)
        
        self.__data = []
        self.__file = ''
        self.__n_data = []
        self.__n_file = ''
        self.__mask_data = []
        self.__mask_file = ''
        
    
    def _read_dna(self, l, lowercase=False):
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
        
        file = os.path.join(self.dir, l.chr + ".dna.2bit")
        
        if not os.path.exists(file):
            return []

        if file != self.__file:
            print('Caching {}...'.format(file))
            self.__file = file
            # Load file into memory
            f = open(file, 'rb')
            self.__data = f.read()
            f.close()
            
            
        return DNA2Bit._read2bit(self.__data, l, offset=True)
    
    
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
        
        d = DNA2Bit._read1bit(self.__n_data, l, offset=True)
        
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
        
        d = DNA2Bit._read1bit(self.__mask_data, l, offset=True)
        
        if mask.startswith(l):
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = DNA_UC_TO_LC_MAP[ret[i]] #ret[i].lower()
        else:
            # Use N as mask
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = DNA_N_UC #'N'
