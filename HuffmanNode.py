# g-sqz
# The node of a Huffman tree

class HuffmanNode:
    # class constructor
    def __init__(self, left, right, data, freq):
        # These values are created
        # when the class is instantiated.
        if (data is None and freq is None):
            self.data = "" + left.data + right.data
            self.freq = 0
            self.left = left
            self.right = right
            self.is_leaf = False
        elif (left is None and right is None):
            self.data = data
            self.freq = freq
            self.left = None
            self.right = None
            self.is_leaf = True

    # basic getters
    def is_leaf(self):
        return self.is_leaf

    def get_left(self):
        return self.left

    def get_right(self):
        return self.right

    def get_data(self):
        return self.data        
