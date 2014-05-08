# g-sqz
# The Huffman tree builder class


import _pickle as pickle
from heapq import *
from struct import *
from HuffmanNode import *


# reads the file and builds a dictionary of data and its frequency
def build_map(file_name):
    if file_name.endswith('.fastq'):
        huffman_map = build_map_fastq(file_name)
    return huffman_map

# builds a map for FASTQ files
def build_map_fastq(file_name):
    huffman_map = {}
    file = open(file_name, 'r')
    read = True
    while read:
        # sequence identifier
        line_1 = file.readline()
        # checks for an empty file or for the end of a file
        if not line_1:
            read = False
        elif line_1[0] != '@':
            raise FileFormatIncorrectException('The sequence identifier line does not match the FASTQ format')
        else:
            # raw sequence
            line_2 = file.readline().rstrip('\n')
            # line 3 is relatively unimportant given it is a comment (for now)
            file.readline()
            # quality scores
            line_4 = file.readline().rstrip('\n')
            # more of an length error check
            if len(line_2) == len(line_4):
                # count the indices in the line (NOTE: does NOT count '\n')
                for i in range(0, len(line_2)):
                    # store the key as sequence-score
                    key = "" + line_2[i] + line_4[i]
                    if key in huffman_map:
                        huffman_map[key] += 1
                    else:
                        huffman_map[key] = 1
            else:
                raise FileFormatIncorrectException('The length of the raw sequence does not match the length of the quality score')
        file.flush()
    file.close()
    return huffman_map

# converts the Huffman map to Huffman nodes
# then builds the Huffman tree
def build_huffman_tree(huffman_map):
    huffman_node_heap = []
    for i, j in huffman_map.items():
        node = HuffmanNode(None, None, i, j)
        heappush(huffman_node_heap, node)
    while (len(huffman_node_heap) > 1):
        left_node = heappop(huffman_node_heap)
        right_node = heappop(huffman_node_heap)
        parent_node = HuffmanNode(left_node, right_node, None, 0)
        heappush(huffman_node_heap, parent_node)
    return heappop(huffman_node_heap)

# generates a dict with <key, val> = <seq-score, huffman_code>
# the dict will provide quick access times while encoding file
def generate_huffman_code_map(huffman_node):
    huffman_code_map = {}
    generate_huffman_code(huffman_node, '', huffman_code_map)
    return huffman_code_map

# recursive function that finds the huffman code
def generate_huffman_code(node, code, huffman_code_map):
    # base case: leaf
    if node.is_leaf():
        huffman_code_map[node.data] = code
    # otherwise: parent node
    else:
        generate_huffman_code(node.left, code+'0', huffman_code_map)
        generate_huffman_code(node.right, code+'1', huffman_code_map)

# generates an exception for invalid file formats
class FileFormatIncorrectException(Exception):
    def __init__(self, error):
        self.error = error
        Exception.__init__(self, 'File Format Incorrect Exception: %s' % error)

# autotest data

# 78 elements
a1 = build_map_fastq('test1.fastq')
b1 = build_huffman_tree(a1)
c1 = generate_huffman_code_map(b1)

# 12159 elements
a2 = build_map_fastq('test2.fastq')
b2 = build_huffman_tree(a2)
c2 = generate_huffman_code_map(b2)

# TODO everything below this <>

def writeHuffFile(filename, newfile):
    hList, tree = buildHuffList(filename)
    # open both files
    new = open(newfile, 'wb')
    old = open(filename, 'r')
    # dump object into new file
    pickle.dump(tree, new)
    nLine = pack('c', '\n')
    new.write(nLine)
    # write the remaining bits
    while True:
        o1 = old.readline()
        if not o1:
            break
        else:
            o2 = old.readline()
            # line 3 is relatively unimportant
            old.readline()
            o4 = old.readline()
            # omits "length = 0" lines ---> pointless!
            if (len(o2) != 1 and len(o4) != 1):
                # read the designated string value
                stt = o1.find(" ") + 1
                stp = o1.find(" ", stt)
                # write the designated string value
                desg = pack('s', o1[stt:stp])
                new.write(desg)
                # write the bits for Huffman
                binCode = ""  # the bitstring
                for i in range(0, len(o2) - 1):
                    nuc = "" + o2[i]
                    val = "" + o4[i]
                    binCode += getBinValue(nuc, val, hList)
                # get and write remainder to file
                rmn = len(binCode) % 8
                rem = chr(rmn)
                new.write(rem)
                # this is epic - makes life a lot easier if DONE HERE
                for i in range(rmn, 8):
                    binCode += "0"
                # intialize x 
                x = 0
                # extracts fragments of 8bits
                for i in range(0, len(binCode)):
                    # extracts indices from fragments
                    c = int(binCode[i])
                    # convert to asci
                    x = x << 1 or c
                    # write bit
                    if (i % 8 == 7):
                        bit = chr(x)
                        new.write(bit)
                        x = 0
                # write new line
                new.write(nLine)
    # close files
    new.close()
    old.close()
    return True


# gets the binary value from the
def getBinValue(n, v, h):
    temp = ['G', 'C', 'A', 'T', 'N']
    subList = h[temp.index(n)]
    for i in range(0, len(subList)):
        if (subList[i][0] == v):
            return subList[i][1]

            # write format:
            # 1. pickled object
            # 2. \n
            # 3. RECURSIVE:
            # a. serial number
            # b. binary remainder
            # c. binary code
            # d. \n

