import sys
import math
import re
import gzip

TWO_BIT_CHAR_MAP = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "a": 0,
    "c": 1,
    "g": 2,
    "t": 3,
    "N": 0,
    "n": 0,
}

FOUR_BIT_CHAR_MAP = {
    "A": 1,
    "C": 2,
    "G": 3,
    "T": 4,
    "a": 5,
    "c": 6,
    "g": 7,
    "t": 8,
    "N": 9,
    "n": 10,
}

def encode_dna2bit(file):
    print(file, file=sys.stderr)
    matcher = re.match(r"(chr(\d+|[XYM]))", file)

    chr = matcher.group(1)

    dna_out = chr + ".dna.2bit"
    mask_out = chr + ".n.1bit"
    repeat_out = chr + ".mask.1bit"

    print("Creating dna files", dna_out, mask_out, repeat_out, "...")

    print("Reading from", file, "...")

    # Extract the sequence as a single line of bases
    if "gz" in file:
        f = gzip.open(file, "rt")
    else:
        f = open(file, "r")

    # skip fasta header
    f.readline()

    sequence = f.readline().strip()
    f.close()

    print("Finished.")

    print("Writing", dna_out, "...")
    fout = open(dna_out, "wb")
    # fout.write(42)

    # How many bytes we need to encode the sequence. We can store
    # 4 bases per byte
    byte_count = int(len(sequence) / 4) + 1

    print("bytes " + str(len(sequence)) + " " + str(byte_count))

    bytes = bytearray(byte_count)

    for i in range(0, len(sequence)):
        # which byte we are in
        bi = int(i / 4)

        # which quarter of the byte to use
        # the first char we encounter must be written to the upper 4 bits
        s = (3 - (i % 4)) * 2

        base = sequence[i]

        if base in TWO_BIT_CHAR_MAP:
            encoded_base = TWO_BIT_CHAR_MAP[base]
        else:
            encoded_base = 0

        # bit shift encoded base into upper or lower half of byte and
        # OR with byte (default value of zero) to encode two bases in one
        # byte
        bytes[bi] = bytes[bi] | (encoded_base << s)

    fout.write(bytes)
    fout.close()

    print("Writing", mask_out, "...")
    fout = open(mask_out, "wb")
    # fout.write(42)

    # How many bytes we need to encode the mask. We can store
    # 8 bases per byte
    byte_count = int(len(sequence) / 8) + 1

    bytes = bytearray(byte_count)

    for i in range(0, len(sequence)):
        # which byte we are in
        bi = int(i / 8)

        # which quarter of the byte to use
        # the first char we encounter must be written to the upper 4 bits
        s = 7 - (i % 8)

        base = sequence[i]

        if base == "N" or base == "n":
            encoded_base = 1
        else:
            encoded_base = 0

        # bit shift encoded base into upper or lower half of byte and
        # OR with byte (default value of zero) to encode two bases in one
        # byte
        bytes[bi] = bytes[bi] | (encoded_base << s)

    fout.write(bytes)
    fout.close()

    print("Writing " + repeat_out + "...\n")
    fout = open(repeat_out, "wb")
    # fout.write(42)

    # How many bytes we need to encode the sequence either upper or lowercase. We can store
    # 8 bases per byte
    byte_count = int(len(sequence) / 8) + 1

    bytes = bytearray(byte_count)

    for i in range(0, len(sequence)):
        # which byte we are in
        bi = int(i / 8)

        # which bit of the byte to use
        # the first char we encounter must be written to the upper 4 bits
        s = 7 - (i % 8)

        base = sequence[i]

        if base == "a" or base == "c" or base == "g" or base == "t" or base == "n":
            encoded_base = 1
        else:
            encoded_base = 0

        # bit shift encoded base into upper or lower half of byte and
        # OR with byte (default value of zero) to encode 8 bases in one
        # byte
        bytes[bi] = bytes[bi] | (encoded_base << s)

    fout.write(bytes)
    fout.close()


def encode_dna4bit(file):
    print(file, file=sys.stderr)
    matcher = re.match(r"(chr(\d+|[XYM]))", file)

    chr = matcher.group(1)

    dna_out = chr.lower() + ".dna.4bit"
    # mask_out = chr + ".n.1bit"
    # repeat_out = chr + ".mask.1bit"

    print("Creating dna files", dna_out, "...")

    print("Reading from", file, "...")

    # Extract the sequence as a single line of bases
    if "gz" in file:
        f = gzip.open(file, "rt")
    else:
        f = open(file, "r")

    # skip fasta header
    f.readline()

    sequence = f.readline().strip()
    f.close()

    print("Finished.")

    print("Writing", dna_out, "...")
    fout = open(dna_out, "wb")
    # first by is 42 for testing endian
    fout.write(bytearray([42]))

    print(f"Sequence is {len(sequence)} bases.")

    # How many bytes we need to encode the sequence. We can store
    # 4 bases per byte
    byte_count = int(len(sequence) / 2) + 1

    print("bytes " + str(len(sequence)) + " " + str(byte_count))

    bytes = bytearray(byte_count)

    for i in range(0, len(sequence)):
        # which byte we are in
        bi = int(i / 2)

        # which half of the byte to use
        # the first char we encounter must be written to the upper 4 bits
        s = (1 - (i % 2)) * 4

        base = sequence[i]

        

        if base in FOUR_BIT_CHAR_MAP:
            encoded_base = FOUR_BIT_CHAR_MAP[base]
        else:
            print("not found", base)
            encoded_base = 0



        # bit shift encoded base into upper or lower half of byte and
        # OR with byte (default value of zero) to encode two bases in one
        # byte
        bytes[bi] = bytes[bi] | (encoded_base << s)

    fout.write(bytes)
    fout.close()
