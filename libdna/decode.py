import os
from abc import ABC, abstractmethod
from typing import Union
import gal

import sys

# from .libdna import parse_loc

# import s3fs

# Use ord('A') etc to get ascii values
DNA_UC_DECODE_DICT = {0: 65, 1: 67, 2: 71, 3: 84}

DNA_LC_DECODE_DICT = {0: 97, 1: 99, 2: 103, 3: 116}
DNA_UC_TO_LC_MAP = {65: 97, 67: 99, 71: 103, 84: 116, 78: 110}

DNA_4BIT_DECODE_DICT = {
    1: 65,
    2: 67,
    3: 71,
    4: 84,
    5: 97,
    6: 99,
    7: 103,
    8: 116,
    9: 78,
    10: 110,
}

DNA_4BIT_COMP_DICT = {
    0: 0,
    65: 84,
    67: 71,
    71: 67,
    84: 65,
    97: 116,
    99: 103,
    103: 99,
    116: 97,
    78: 78,
    110: 10,
}


DNA_N_UC = 78
DNA_N_LC = 110
DNA_COMP_DICT = {
    65: 84,
    67: 71,
    84: 65,
    71: 67,
    97: 116,
    99: 103,
    116: 97,
    103: 99,
    78: 78,
    110: 110,
}
EMPTY_BYTEARRAY = bytearray(0)

SHIFT_4BIT_MAP = {0: 4, 1: 0}
SHIFT_2BIT_MAP = {0: 6, 1: 4, 2: 2, 3: 0}
SHIFT_1BIT_MAP = {0: 7, 1: 6, 2: 5, 3: 4, 4: 3, 5: 2, 6: 1, 7: 0}


class DNA(ABC):
    @abstractmethod
    def dna(self, *args):
        raise NotImplementedError

    def fasta(self, loc: gal.genomic.Location, mask="upper"):
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

        # l = libdna.parse_loc(loc)

        print(f">{loc}")
        print(self.dna(loc, mask=mask))

    def to_fasta(self, loc: gal.genomic.Location, mask="upper"):
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

        # l = libdna.parse_loc(loc)

        return f">{loc}\n{self.dna(loc, mask=mask)}"

    def to_file(self, loc: gal.genomic.Location, file, mask: str = "upper"):
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

        f = open(file, "w")
        print(f">{loc}", file=f)
        print(f"{self.dna(loc, mask=mask)}", file=f)
        f.close()


class DNAStr(DNA):
    def __init__(self, dir):
        self._dir = dir

    def dna(self, loc: gal.genomic.Location):
        # l = libdna.parse_loc(loc)

        sstart = loc.start - 1
        send = loc.end - 1

        file = os.path.join(self._dir, loc.chr + ".txt")

        f = open(file, "r")

        f.seek(sstart)

        length = send - sstart + 1

        seq = f.read(length)

        f.close()

        return seq


class DNABin(DNA):
    def read_data(self, file: str, seek: int, n: int) -> bytes:
        """
        Reads data from a file source

        Parameter
        ---------
        file : str
            Relative path to file
        seek : int
            Start offset in bytes
        n : int
            Amount of data to read in bytes

        Returns
        -------
        bytearray
            Data from file
        """
        file = os.path.join(self.dir, file).lower()

        if not os.path.exists(file):
            return None

        #print(seek)

        with open(file, "rb") as f:
            # first byte is 42
            f.seek(seek)
            data = f.read(n)  # l.length // 8 + 2)

        return data


class DNA2Bit(DNABin):
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

        l = int(len(dna) / 2)

        for i in range(0, l):
            b = DNA_COMP_DICT[dna[i]]
            dna[i] = DNA_COMP_DICT[dna[i2]]
            dna[i2] = b
            i2 -= 1

    def _read1bit(self, d: bytes, loc: gal.genomic.Location, offset=False) -> bytearray:
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

        ret = bytearray([0] * loc.length)

        if offset:
            bi = int(s / 8)
        else:
            bi = 0

        for i in range(0, length):
            block = s % 8

            v = d[bi] >> SHIFT_1BIT_MAP[block]

            if block == 7:
                bi += 1

            # if block == 0:
            #     v = (d[bi] >> 7)
            # elif block == 1:
            #     v = (d[bi] >> 6)
            # elif block == 2:
            #     v = (d[bi] >> 5)
            # elif block == 3:
            #     v = (d[bi] >> 4)
            # elif block == 4:
            #     v = (d[bi] >> 3)
            # elif block == 5:
            #     v = (d[bi] >> 2)
            # elif block == 6:
            #     v = (d[bi] >> 1)
            # else:
            #     v = d[bi]
            #     bi += 1

            # Only care about the lowest bit
            ret[i] = v & 1

            s += 1

        return ret

    def _read2bit(self, d: bytes, loc: gal.genomic.Location, offset=False) -> bytearray:
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

        # print(d, loc)

        s = loc.start - 1

        ret = bytearray([0] * loc.length)  # []

        if offset:
            bi = int(s / 4)
        else:
            bi = 0

        for i in range(loc.length):
            block = s % 4

            v = d[bi] >> SHIFT_2BIT_MAP[block]

            if block == 3:
                bi += 1

            # if block == 0:
            #     v = (d[bi] >> 6)
            # elif block == 1:
            #     v = (d[bi] >> 4)
            # elif block == 2:
            #     v = (d[bi] >> 2)
            # else:
            #     v = d[bi]

            #     # Reached end of byte so we are moving into the next byte
            #     bi += 1

            # Only care about the lowest 2 bits (3)
            ret[i] = DNA_UC_DECODE_DICT[v & 3]

            s += 1

        return ret

    def _read_dna(self, loc: gal.genomic.Location, lowercase=False) -> bytearray:
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

        file = f"{loc.chr}.dna.2bit"

        # print(loc.length, file=sys.stderr)

        s = loc.start - 1
        e = s + loc.length
        bs = int(s / 4)
        be = int(e / 4)
        l = be - bs + 1

        data = self.read_data(file, bs, l)

        # print(loc, data, file=sys.stderr)

        return self._read2bit(data, loc)

    def _read_1bit_file(self, file: str, loc: gal.genomic.Location):
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
        bs = int(s / 8)
        be = int(e / 8)
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

    def _read_n(self, loc: gal.genomic.Location, ret: bytearray):
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

        file = f"{loc.chr}.n.1bit"

        data = self._read_1bit_file(file, loc)

        d = self._read1bit(data, loc)

        for i in range(0, len(ret)):
            if d[i] == 1:
                ret[i] = DNA_N_UC  # 'N'

    def _read_mask(self, loc: gal.genomic.Location, ret, mask="upper"):
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

        if mask.startswith("u"):
            return

        file = "{}.mask.1bit".format(loc.chr)

        data = self._read_1bit_file(file, loc)

        d = self._read1bit(data, loc)

        if mask.startswith("l"):
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = DNA_UC_TO_LC_MAP[ret[i]]  # ret[i].lower()
        else:
            # Use N as mask
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = DNA_N_UC  # 'N'

    def dna(
        self, loc: gal.genomic.Location, mask="lower", rev_comp=False, lowercase=False
    ):
        """
        Returns the DNA for a location.

        Parameters
        ----------
        loc : libdna.Loc
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

        ret = self._read_dna(loc, lowercase=lowercase)

        self._read_n(loc, ret)

        self._read_mask(loc, ret, mask=mask)

        if rev_comp:
            DNA2Bit._rev_comp(ret)

        ret = ret.decode("utf-8")

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

        s1 = r1.pos  # + 1

        # end of first read
        e1 = s1 + r1.length - 1

        # start of second read
        s2 = r2.pos  # + 1

        e2 = s2 + r2.length - 1

        inner = s2 - e1 - 1

        if inner >= 0:
            # get the dna sequence spanning the two reads
            seq = self.dna(gal.genomic.Location(r1.chr, s1, e2))
        else:
            # Reads overlap so concatenate the first read with the
            # portion of the second read that is not overlapping
            # (inner is negative so flip sign for array indexing)
            seq = r1.seq + r2.seq[-inner:]

        return seq


class DNA4Bit(DNABin):
    def __init__(self, dir):
        self._dir = dir

    @property
    def dir(self):
        return self._dir

    @staticmethod
    def rev_comp(dna):
        """
        Parameters
        ----------
        dna : bytearray
            dna sequence to be reverse complemented
        """

        i2 = len(dna) - 1

        l = int(len(dna) / 2)

        for i in range(0, l):
            b = DNA_4BIT_COMP_DICT[dna[i]]
            dna[i] = DNA_4BIT_COMP_DICT[dna[i2]]
            dna[i2] = b
            i2 -= 1

    def _read4bit(self, d: bytes, loc: gal.genomic.Location, offset=False) -> bytearray:
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

        # print(d, loc)

        s = loc.start - 1

        ret = bytearray([0] * loc.length)  # []

        if offset:
            bi = int(s / 2)
        else:
            bi = 0

        #print([d[bi] for bi in range(len(d))])

        for i in range(loc.length):
            block = s % 2

            v = d[bi] >> SHIFT_4BIT_MAP[block]

            #print(i, bi, d[bi], v)

            # if block == 0:
            #     v = (d[bi] >> 6)
            # elif block == 1:
            #     v = (d[bi] >> 4)
            # elif block == 2:
            #     v = (d[bi] >> 2)
            # else:
            #     v = d[bi]

            #     # Reached end of byte so we are moving into the next byte
            #     bi += 1

            # Only care about the lowest 4 bits (3)
            ret[i] = DNA_4BIT_DECODE_DICT[v & 15]

            if block == 1:
                bi += 1

            s += 1

        return ret

    def _read_dna(self, loc: gal.genomic.Location, lowercase=False) -> bytearray:
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

        file = f"{loc.chr}.dna.4bit"

        # print(loc.length, file=sys.stderr)

        s = loc.start - 1
        e = s + loc.length
        bs = int(s / 2)
        be = int(e / 2)
        l = be - bs + 1

        # skip first byte as this is 42
        data = self.read_data(file, bs + 1, l)
     
        return self._read4bit(data, loc)

    def dna(
        self, loc: gal.genomic.Location, mask="lower", rev_comp=False, lowercase=False
    ):
        """
        Returns the DNA for a location.

        Parameters
        ----------
        loc : libdna.Loc
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

        ret = self._read_dna(loc, lowercase=lowercase)

        if rev_comp:
            DNA4Bit._rev_comp(ret)

        ret = ret.decode("utf-8")

        if lowercase:
            ret = ret.lower()

        return ret


# class S3DNA2Bit(DNA2Bit):
#     def __init__(self, bucket, dir):
#         super().__init__(dir)
#         self.__bucket = bucket
#         self.__fs = s3fs.S3FileSystem(anon=True)

#     @property
#     def bucket(self):
#         return self.__bucket

#     def read_data(self, file, seek, n):
#         file = os.path.join(self.__bucket, self.dir, file).lower()

#         f = self.__fs.open(file, 'rb')
#         f.seek(seek)
#         data = f.read(n)
#         return data


class CachedDNA2Bit(DNA2Bit):
    def __init__(self, dir):
        super().__init__(dir)

        self.__data = []
        self.__file = ""
        self.__n_data = []
        self.__n_file = ""
        self.__mask_data = []
        self.__mask_file = ""

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
            print(file, "does not exist.")
            return EMPTY_BYTEARRAY

        if file != self.__file:
            print("Caching {}...".format(file))
            self.__file = file
            # Load file into memory
            f = open(file, "rb")
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
            print("Caching {}...".format(file))
            f = open(file, "rb")
            self.__n_data = f.read()
            f.close()
            self.__n_file = file

        d = self._read1bit(self.__n_data, l, offset=True)

        for i in range(0, len(ret)):
            if d[i] == 1:
                ret[i] = DNA_N_UC  # 'N'

    def _read_mask(self, l, ret, mask="upper"):
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

        if mask.startswith("u"):
            return

        file = os.path.join(self.dir, l.chr + ".mask.1bit")

        if not os.path.exists(file):
            return

        if file != self.__mask_file:
            print("Caching {}...".format(file))
            f = open(file, "rb")
            self.__mask_data = f.read()
            f.close()
            self.__mask_file = file

        d = self._read1bit(self.__mask_data, l, offset=True)

        if mask.startswith("l"):
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = DNA_UC_TO_LC_MAP[ret[i]]  # ret[i].lower()
        else:
            # Use N as mask
            for i in range(0, len(ret)):
                if d[i] == 1:
                    ret[i] = DNA_N_UC  # 'N'
