# g-sqz
# The Huffman tree builder class

from heapq import *
from _pickle import *
from HuffmanNode import *

# reads the file and builds a dictionary of data and its frequency
def build_map(file_name):
    if file_name.endswith('.fastq'):
        huffman_map = build_map_fastq(file_name)
    return huffman_map

# builds a map for FASTQ files
def build_map_fastq(file_name):
    huffman_map = {}
    line_1_list = []
    file = open(file_name, 'r')
    line_1 = file.readline()
    line_len = 0
    while line_1:
        if line_1[0] != '@':
            raise FileFormatIncorrectException('The sequence identifier line does not match the FASTQ format')
        else:
            line_1_list.append(line_1)
            # raw sequence
            line_2 = file.readline().rstrip('\n')
            # line 3 - TBD
            file.readline()
            # quality scores
            line_4 = file.readline().rstrip('\n')
            # more of an length error check
            if len(line_2) == len(line_4):
                if line_len == 0:
                    line_len = len(line_2)
                # count the indices in the line (NOTE: does NOT count '\n')
                for i in range(len(line_2)):
                    # store the key as sequence-score
                    key = '' + line_2[i] + line_4[i]
                    if key in huffman_map:
                        huffman_map[key] += 1
                    else:
                        huffman_map[key] = 1
            else:
                raise FileFormatIncorrectException('The length of the raw sequence does not match the length of the quality score')
            file.flush()
            line_1 = file.readline()
    file.close()
    return huffman_map, line_1_list, line_len


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
    huffman_map, line_1_list, line_len = build_map(file_name)
    huffman_node = build_huffman_tree(huffman_map)
    huffman_encode_map = generate_huffman_code_map(huffman_node)
    huffman_decode_map = {val:key for key, val in huffman_encode_map.items()}
    
    # file io
    read_file = open(file_name, 'r')
    gsqz_name = file_name+'.gsqz'
    write_file = open(gsqz_name, 'wb')
    
    # 1. write line length
    write_file.write(line_len.to_bytes(1, byteorder='big'))
    
    # 2. dump objects and their lengths
    pickled_list = dumps(line_1_list)
    write_file.write(len(pickled_list).to_bytes(3, byteorder='big'))
    pickled_map = dumps(huffman_decode_map)
    write_file.write(len(pickled_map).to_bytes(3, byteorder='big'))
    write_file.write(pickled_list)
    write_file.write(pickled_map)
    write_file.close()
    
    # 3. write lines
    raw_code = ''
    line_1 = read_file.readline()
    while line_1:
        line_2 = read_file.readline().rstrip('\n')
        read_file.readline()
        line_4 = read_file.readline().rstrip('\n')            
        for i in range(len(line_2)):
            seq_scr = '' + line_2[i] + line_4[i]
            raw_code += huffman_encode_map[seq_scr]
        rem = -(len(raw_code)%8)
        append_bytes(gsqz_name, raw_code[:rem])
        raw_code = raw_code[rem:]
        read_file.flush()
        line_1 = read_file.readline()
    read_file.close()
    rem = 8-len(raw_code)
    if rem > 0:
        raw_code += '0'*rem
        append_bytes(gsqz_name, raw_code)
    print('Successfully encoded: ' + file_name)
    return huffman_map, line_1_list, line_len, huffman_node, huffman_encode_map

# appends bytes to output file
def append_bytes(file_name, bin_str):    
    byte_bin_map = byte_bin(False)
    byte_str = b''
    write_file = open(file_name, 'ab')
    for i in range(0, len(bin_str), 8):
        #print(bin_str[i:i+8])
        byte_str += byte_bin_map[bin_str[i:i+8]]
    write_file.write(byte_str)
    write_file.close()    
    
# byte-binary map
def byte_bin(bytetobin=True):
    byte_bin_map = {}
    if bytetobin:
        for i in range(256):
            bin_val = bin(i)[2:]
            prefix_bin_val = 8-len(bin_val)
            if (prefix_bin_val > 0):
                bin_val = '0'* prefix_bin_val + bin_val
            byte_val = i.to_bytes(1, byteorder='big')
            byte_bin_map[byte_val] = bin_val
    else:
        for i in range(256):
            bin_val = bin(i)[2:]
            prefix_bin_val = 8-len(bin_val)
            if (prefix_bin_val > 0):
                bin_val = '0'* prefix_bin_val + bin_val
            byte_val = i.to_bytes(1, byteorder='big')
            byte_bin_map[bin_val] = byte_val
    return byte_bin_map

# TODO: complete function <high>
# decodes gsqz file        
def gsqz_decode_fastq_simple(file_name):    
    byte_bin_map = byte_bin()
    write_file = open(file_name+'.fastq', 'w')
    read_file = open(file_name, 'rb')
    pickled_map_raw_len = ''
    for i in range(3):
        pickled_map_raw_len += byte_bin_map[read_file.read(1)]
    pickled_map_len = int(pickled_map_raw_len, 2)
    print(pickled_map_len)        
    pickled_map = loads(read_file.read(pickled_map_len))
    print(pickled_map)
    max_line_len = read_file.read(1)
    curr_len = 0
    char_str = ''
    seq = ''
    scr = ''
    byte = read_file.read(1)
    while byte:
        char_str += byte_bin_map[byte]
        byte = read_file.read(1)
    found = False
    while not found:
        char_str_len = len(char_str)
        for i in range(char_str_len):
            for j in range(1, char_str_len-i):
                key = char_str[i:i+j]
                if key in byte_bin_map:
                    seq_scr = byte_bin_map[key]
                    seq += seq_scr[0]
                    scr += seq_scr[1]
    
# generates an exception for invalid file formats
class FileFormatIncorrectException(Exception):
    def __init__(self, error):
        self.error = error
        Exception.__init__(self, 'File Format Incorrect Exception: %s' % error)

# autotest data
if __name__ == '__main__':
    # 78 elements
    a1, a2, a3, a4, a5 = gsqz_encode_fastq_simple('test1.fastq')

    # 12159 elements
    b1, b2, b3, b4, b5 = gsqz_encode_fastq_simple('test2.fastq')
