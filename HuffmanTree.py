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
                    key = '' + line_2[i] + line_4[i]
                    if key in huffman_map:
                        huffman_map[key] += 1
                    else:
                        huffman_map[key] = 1
            else:
                raise FileFormatIncorrectException('The length of the raw sequence does not match the length of the quality score')
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

# writes the g-sqz'd file
def gsqz_encode_fastq_simple(file_name):
    # prepare gsqz data
    huffman_map = build_map(file_name)
    huffman_node = build_huffman_tree(huffman_map)
    huffman_encode_map = generate_huffman_code_map(huffman_node)
    huffman_decode_map = {val:key for key, val in huffman_encode_map.items()}
    # write data to file
    read_file = open(file_name, 'r')
    write_file = open(file_name+'.gsqz', 'wb')
    # 1. dump pickled map
    pickle.dump(huffman_decode_map, write_file)
    read = True
    while read:
        line_1 = read_file.readline()
        if not line_1:
            read = False
        elif line_1[0] != '@':
            raise FileFormatIncorrectException('The sequence identifier line does not match the FASTQ format')
        else:
            line_2 = read_file.readline().rstrip('\n')
            read_file.readline()
            line_4 = read_file.readline().rstrip('\n')
            if len(line_2) == len(line_4):
                # generate raw encode string
                raw_code = ''
                for i in range(0, len(line_2)):
                    seq_scr = '' + line_2[i] + line_4[i]
                    raw_code += huffman_encode_map[seq_scr]
                # add 0s to end
                remainder = len(raw_code)%8
                if remainder != 0:
                    balance = 8 - remainder
                    for i in range(balance):
                        raw_code += '0'
                print(raw_code)
                # write corresponding bytes to file
                byte_str = ''
                for i in range(0, len(raw_code), 8):
                    int_val = int(raw_code[i:i+8], 2)
                    write_file.write(int_val.to_bytes(1, byteorder='big'))                                     
            else:
                raise FileFormatIncorrectException('The length of the raw sequence does not match the length of the quality score')
        write_file.flush()
    read_file.close()
    write_file.close()

# generates an exception for invalid file formats
class FileFormatIncorrectException(Exception):
    def __init__(self, error):
        self.error = error
        Exception.__init__(self, 'File Format Incorrect Exception: %s' % error)

# autotest data

# 78 elements
gsqz_encode_fastq_simple('test1.fastq')
a1 = build_map_fastq('test1.fastq')
b1 = build_huffman_tree(a1)
c1 = generate_huffman_code_map(b1)

# 12159 elements
#gsqz_encode_fastq_simple('test2.fastq')
a2 = build_map_fastq('test2.fastq')
b2 = build_huffman_tree(a2)
c2 = generate_huffman_code_map(b2)
