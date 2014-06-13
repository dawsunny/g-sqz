# g-sqz
# The Huffman tree builder class

from collections import deque
from heapq import *
from _pickle import *
from HuffmanNode import *
import gc
import os
import re

# reads the file and builds a dictionary of data and its frequency
def build_map(file_name):
    if file_name.endswith('.fastq'):
        huffman_map = build_map_fastq(file_name)
    return huffman_map

# builds a map for FASTQ files
def build_map_fastq(file_name):
    huffman_map = {}
    seq_pattern = []
    separators = '@_.:-'
    file = open(file_name, 'r')
    line_1 = file.readline().lstrip('@').rstrip('\n')
    line_len = 0

    sep_list = []
    str_list = []
    list_pos = None
    if line_1:
        if ' length' in line_1:
            line_1 = line_1[:line_1.index(' length')]
        sep_list = re.findall(r'[ _/\.,:;#~-]+', line_1)
        str_list = re.split(r'[ _/\.,:;#~-]+', line_1)
        list_pos = [None]*len(str_list)
    while line_1:
        str_list_temp = re.split(r'[ _/\.,:;#~-]+', line_1)
        for i in range(len(str_list)):
            str_temp = str_list[i]
            if str_temp != '' and str_temp != str_list_temp[i]:
                max_pos = str_max_pos(str_temp, str_list_temp[i])
                if max_pos < len(str_temp):
                    str_list[i] = str_list[i][:max_pos]
                    list_pos[i] = max_pos
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
        line_1 = file.readline().lstrip('@').rstrip('\n')
    file.close()
    str_opt, str_pos = optimize_seq(sep_list, str_list, list_pos)    
##    print(sep_list, str_list, str_opt)
##    print(list_pos, str_pos)
    return huffman_map, line_len, str_opt, str_pos, list_pos


# returns the maximum length of the substring that matches
def str_max_pos(string_1, string_2):
    start_val = min(len(string_1), len(string_2))
    max_pos = 0
    for i in range(start_val, 0, -1):
        if string_1[:i] == string_2[:i]:
            max_pos = i
            break
    return max_pos


# returns the string made of the common characters in the SEQ line
def optimize_seq(sep_list, str_list, list_pos):
    str_count = 0
    str_pos = []
    str_opt = ''
    sep_len = len(sep_list)
    for i in range(len(list_pos)):
        pos_temp = list_pos[i]
        str_temp = str_list[i]
        if pos_temp is None:
            str_opt += str_temp
            str_count += len(str_list[i])
        elif pos_temp == 0:
            str_pos.append(str_count)
        else:
            str_opt += str_temp
            str_count += len(str_list[i])
            str_pos.append(str_count)
        if i < sep_len:
            str_opt += sep_list[i]
            str_count += 1
    return str_opt, str_pos


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
def gsqz_encode_fastq(file_name):
    # preparation
    huffman_map, line_len, str_opt, str_pos, list_pos = build_map(file_name)
    huffman_node = build_huffman_tree(huffman_map)
    huffman_encode_map = generate_huffman_code_map(huffman_node)
    huffman_decode_map = {val:key for key, val in huffman_encode_map.items()}
    huffman_decode_map[-1] = (str_opt, tuple(str_pos), tuple(list_pos))
    seek_map = {}
    
    # file io
    read_file = open(file_name, 'r')
    gsqz_name = file_name+'.gsqz'
    temp_name = file_name+'.tmp'
    write_file = open(gsqz_name, 'wb')
    
    line_1 = read_file.readline().lstrip('@').rstrip('\n')
    line_end = len(line_1)
    if line_1:
        if ' length' in line_1:
            line_len += 128
            line_end = line_1.index(' length')
    
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
    
    while line_1:
        if ' length' in line_1:
            line_1 = line_1[:line_end]
        str_list_temp = re.split(r'[ _/\.,:;#~-]+', line_1)
        line_1_keys = []
        for i in range(len(list_pos)):
            pos_temp = list_pos[i]
            if pos_temp == 0:
                line_1_keys.append(int(str_list_temp[i]))
            elif pos_temp is not None:
                line_1_keys.append(int(str_list_temp[i][pos_temp:]))
        seek_map[tuple(line_1_keys)] = (byte_index, bit_index)
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
        line_1 = read_file.readline().lstrip('@').rstrip('\n')
    read_file.close()
    if len(raw_code) > 0:
        raw_code += '0'*(8-len(raw_code))
        append_bytes(temp_name, raw_code)

    # dump seek map
    write_file = open(gsqz_name, 'ab')
    pickled_seek = dumps(seek_map)
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
    fastq_size = os.path.getsize(file_name)/1024
    gsqz_size = os.path.getsize(gsqz_name)/1024
    print(file_name)
    print('Seek dict size: {:.2f}KB'.format(len(pickled_seek)/1024))
    print('Decode dict size: {:.2f}KB'.format(len(pickled_decode)/1024))
    print('Uncompressed Size: {:.2f}KB'.format(fastq_size))
    print('Compressed Size: {:.2f}KB'.format(gsqz_size))
    print('Compression Savings: {:.2f}'.format(1-(gsqz_size/fastq_size)))
    print('---')
    return huffman_map, huffman_node, huffman_encode_map, line_len, seek_map, str_opt, list_pos, str_pos

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
def gsqz_decode_fastq(gsqz_file, decode='full', start=None, end=None):
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

    # prepare seek map
    min_huffman = len(min(huffman_decode_map.keys(), key=len))
    
    stt_val = None
    end_val = None
    if decode == 'full':
        stt_val = (0, 0)
        end_val = (float('inf'), 0)          
    elif start in seek and end in seek:
        stt_val = seek[start]
        end_val = seek[end]
    else:
        raise KeyError()

    seq_heap = []
    for key,val in seek.items():
        if stt_val <= val:
            heappush(seq_heap, (val[0], key))

    # read and convert entire file to a huffman code string
    char_str = ''
    if stt_val != (0, 0):
        read_file.seek(stt_val[0]-1, 1)
        byte = read_file.read(1)
        char_str += byte_bin_map[byte][stt_val[1]:]        
    
    byte = read_file.read(1)
    while byte:
        char_str += byte_bin_map[byte]
        byte = read_file.read(1)
    read_file.close()

    # decode huffman code string
    if decode == 'full':        
        stt_pos = 0
        end_pos = min_huffman
        while len(seq_heap) > 0:
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
            append_block(write_file, heappop(seq_heap)[1], seq, scr)

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
def main():
    rel = 5
    ver = 3

    for i in range(1, rel):
        for j in range(ver):
            file_name = './test/test{:d}.{:d}.fastq'.format(i, j)
            gsqz_encode_fastq(file_name)
            gc.collect()
            #gsqz_decode_fastq(file_name)
            gc.collect()

    #gsqz_encode_fastq('./test/test3.2.fastq')

if __name__ == '__main__':
    main()
