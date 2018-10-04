import sys
import math
import re

CHAR_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3, 'N': 0, 'n': 0}


def encode_dna2bit(file):
    matcher = re.match(r'(chr(\d+|[XYM]))', file)

    chr = matcher.group(1)

    dna_out = chr + "dna.2bit"
    mask_out = chr + ".n.1bit"
    repeat_out = chr + ".mask.1bit"

    print("Creating dna files", dna_out, mask_out, repeat_out, "...")


    print("Reading from " + file + "...\n")

    # Extract the sequence as a single line of bases
    f = open(file, 'r')
    sequence = f.readline().strip()
    f.close()

    print("Finished.\n")



    print("Writing", dna_out, "...")
    fout = open(dna_out, 'wb')

    # How many bytes we need to encode the sequence. We can store
    # 4 bases per byte
    byte_count = int(math.ceil(len(sequence) / 4))

    sys.stderr.write("byte " + str(len(sequence)) + " " + str(byte_count) + " \n")

    bytes = bytearray(byte_count)

    for i in range(0, len(sequence)):
        # which byte we are in    
        bi = i // 4
      
        # which quarter of the byte to use
        # the first char we encounter must be written to the upper 4 bits
        s = (3 - (i % 4)) * 2

        base = sequence[i]
      
        if base in CHAR_MAP:
            encoded_base = CHAR_MAP[base]
        else:
            encoded_base = 0
      
        # bit shift encoded base into upper or lower half of byte and
        # OR with byte (default value of zero) to encode two bases in one
        # byte
        bytes[bi] = bytes[bi] | (encoded_base << s)

    fout.write(bytes)
    fout.close()


    print("Writing", mask_out, "...")
    fout = open(mask_out, 'wb')

    # How many bytes we need to encode the sequence. We can store
    # 4 bases per byte
    byte_count = int(math.ceil(len(sequence) / 8))

    bytes = bytearray(byte_count)

    for i in range(0, len(sequence)):
        # which byte we are in    
        bi = i // 8
      
        # which quarter of the byte to use
        # the first char we encounter must be written to the upper 4 bits
        s = (7 - (i % 8))

        base = sequence[i]
      
        if base == 'N' or base == 'n':
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
    fout = open(repeat_out, 'wb')

    # How many bytes we need to encode the sequence. We can store
    # 4 bases per byte
    byte_count = int(math.ceil(len(sequence) / 8))

    bytes = bytearray(byte_count)

    for i in range(0, len(sequence)):
        # which byte we are in    
        bi = i // 8
      
        # which quarter of the byte to use
        # the first char we encounter must be written to the upper 4 bits
        s = (7 - (i % 8))

        base = sequence[i]
      
        if base == 'a' or base == 'c' or base == 'g' or base == 't' or base == 'n':
            encoded_base = 1
        else:
            encoded_base = 0
      
        # bit shift encoded base into upper or lower half of byte and
        # OR with byte (default value of zero) to encode two bases in one
        # byte
        bytes[bi] = bytes[bi] | (encoded_base << s)

    fout.write(bytes)
    fout.close()
