# g-sqz
# The node of a Huffman tree

class HuffmanNode:
    # class constructor
    # better to have a list of data than a concatenated string in case on DNA represented numerically
    def __init__(self, left, right, data, freq):
        # leaf node
        if (left is None or right is None):
            self.left = None
            self.right = None
            self.data = data
            self.freq = freq
        # parent node
        else:
            self.left = left
            self.right = right
            self.data = None
            self.freq = left.freq + right.freq

    # for comparison (heapq)
    def __lt__(self, other):
        return self.freq < other.freq
    
    def __le__(self, other):
        return self.freq <= other.freq
    
    def __eq__(self, other):
        return self.freq == other.freq
    
    def __ne__(self, other):
        return self.freq != other.freq
    
    def __gt__(self, other):
        return self.freq > other.freq
    
    def __ge__(self, other):
        return self.freq >= other.freq
    
    # basic getters
    def is_leaf(self):
        return (self.left is None and self.right is None)    
