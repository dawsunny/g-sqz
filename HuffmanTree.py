# g-sqz
# The Huffman tree builder class

from heapq import *
from _pickle import *
from HuffmanNode import *
import os

# reads the file and builds a dictionary of data and its frequency
def build_map(file_name):
    if file_name.endswith('.fastq'):
        huffman_map = build_map_fastq(file_name)
    return huffman_map

# builds a map for FASTQ files
def build_map_fastq(file_name):
    huffman_map = {}
    file = open(file_name, 'r')
    line_1 = file.readline()
    line_len = 0
    while line_1:
        if line_1[0] != '@':
            raise FileFormatIncorrectException('The sequence identifier line does not match the FASTQ format')
        else:
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
    return huffman_map, line_len


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
def gsqz_encode_fastq(file_name, encode='basic'):
    # preparation
    huffman_map, line_len = build_map(file_name)
    huffman_node = build_huffman_tree(huffman_map)
    huffman_encode_map = generate_huffman_code_map(huffman_node)
    huffman_decode_map = {val:key for key, val in huffman_encode_map.items()}
    if encode == 'basic':
        seek_v1 = []
    elif encode == 'mapped':
        seek_map = {}
    else:        
        seek_v2 = []   
    
    # file io
    read_file = open(file_name, 'r')
    gsqz_name = file_name+'.gsqz'
    temp_name = file_name+'.tmp'
    write_file = open(gsqz_name, 'wb')
    
    # write line length
    write_file.write(line_len.to_bytes(1, byteorder='big'))
    
    # dump Huffman map
    pickled_decode = dumps(huffman_decode_map)
    write_file.write(len(pickled_decode).to_bytes(3, byteorder='big'))
    write_file.write(pickled_decode)
    write_file.close()
    
    # write bytes to temp file, build seek map
    raw_code = ''
    byte_index = 0
    bit_index = 0
    line_1 = read_file.readline()
    while line_1:
        if encode == 'basic':
            seek_v1.append(line_1.strip('@\n'))
        elif encode == 'mapped':
            seek_map[line_1.strip('@\n')] = (byte_index, bit_index)
        else:        
            seek_v2.append((byte_index, bit_index))
        line_2 = read_file.readline().rstrip('\n')
        read_file.readline()
        line_4 = read_file.readline().rstrip('\n')            
        for i in range(len(line_2)):
            seq_scr = '' + line_2[i] + line_4[i]
            raw_code += huffman_encode_map[seq_scr]
        raw_code_len = len(raw_code)
        rem = raw_code_len % 8
        bit_index = rem
        byte_index += raw_code_len // 8
        if rem == 0:
            append_bytes(temp_name, raw_code)
            raw_code = ''
        else:
            append_bytes(temp_name, raw_code[:-rem])
            raw_code = raw_code[-rem:]
        read_file.flush()
        line_1 = read_file.readline()
    read_file.close()
    if len(raw_code) > 0:
        raw_code += '0'*(8-len(raw_code))
        append_bytes(temp_name, raw_code)

    # dump seek map
    write_file = open(gsqz_name, 'ab')
    seek = None
    if encode == 'basic':
        seek = seek_v1
    elif encode == 'mapped':
        seek = seek_map
    else:        
        seek = seek_v2
    pickled_seek = dumps(seek)
    write_file.write(len(pickled_seek).to_bytes(3, byteorder='big'))
    write_file.write(pickled_seek)
    
    # append temp file bytes, delete temp file
    temp_file = open(temp_name, 'rb')
    write_file.write(temp_file.read())   
    write_file.close()
    temp_file.close()

    # delete temp file
    os.remove(temp_name)
    
    # confirmation and output
    print('Successfully encoded: ' + file_name)
    return huffman_map, huffman_node, huffman_encode_map, line_len

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


# decodes gsqz file        
def gsqz_decode_fastq(gsqz_file, decode='basic'):
    # preparation
    byte_bin_map = byte_bin()
    read_file = open(gsqz_file, 'rb')    
    write_file = gsqz_file + '.fastq'
    
    # delete old file if it exists since decode appends in chunks
    try:
        os.remove(write_file)
    except OSError:
        pass    

    # decode metadata from file
    line_len = int.from_bytes(read_file.read(1), byteorder='big')
    pickled_map_len = int.from_bytes(read_file.read(3), byteorder='big')
    huffman_decode_map = loads(read_file.read(pickled_map_len))
    pickled_seek_len = int.from_bytes(read_file.read(3), byteorder='big')
    seek = loads(read_file.read(pickled_seek_len))
    a = huffman_decode_map

    # read and convert entire file to a huffman code string
    char_str = ''
    byte = read_file.read(1)
    while byte:
        char_str += byte_bin_map[byte]
        byte = read_file.read(1)
    read_file.close()

    # decode huffman code string
    if decode == 'basic':
        min_huffman = len(min(huffman_decode_map.keys(), key=len))
        stt_pos = 0
        end_pos = min_huffman
        for i in seek:
            curr_len = 0
            seq = ''
            scr = ''
            while curr_len < line_len:
                seq_scr_raw = char_str[stt_pos:end_pos]
                if seq_scr_raw in huffman_decode_map:
                    seq_scr = huffman_decode_map[seq_scr_raw]
                    seq += seq_scr[0]
                    scr += seq_scr[1]
                    curr_len += 1
                    stt_pos = end_pos
                    end_pos = stt_pos + min_huffman
                else:
                    end_pos += 1
                    
            # append
            append_block(write_file, i, seq, scr)

    # confirmation
    print('Successfully decoded: ' + gsqz_file)
    return line_len, pickled_map_len, huffman_decode_map, pickled_seek_len, seek
        
# appends block to output file
def append_block(file_name, line_1, line_2, line_4):
    nwl = '\n'
    file = open(file_name, 'a')
    file.write('@' + line_1 + nwl)
    file.write(line_2 + nwl)
    file.write('+' + nwl)
    file.write(line_4 + nwl)
    file.close()
    
# generates an exception for invalid file formats
class FileFormatIncorrectException(Exception):
    def __init__(self, error):
        self.error = error
        Exception.__init__(self, 'File Format Incorrect Exception: %s' % error)

# autotest data
if __name__ == '__main__':
    # 78 elements
    a1, a2, a3, a4 = gsqz_encode_fastq('./test/test0.0.fastq')

    # 12159 elements
    b1, b2, b3, b4 = gsqz_encode_fastq('./test/test0.1.fastq')

    #decode
    c1, c2, c3, c4, c5 = gsqz_decode_fastq('./test/test0.0.fastq.gsqz')
    d1, d2, d3, d4, d5 = gsqz_decode_fastq('./test/test0.1.fastq.gsqz')
