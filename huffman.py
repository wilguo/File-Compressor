"""
Code for compressing and decompressing using Huffman compression.
"""

from nodes import HuffmanNode, ReadNode

# ====================
# Helper functions for manipulating bytes


def get_bit(byte, bit_num):
    """ Return bit number bit_num from right in byte.
    @param int byte: a given byte
    @param int bit_num: a specific bit number within the byte
    @rtype: int
    >>> get_bit(0b00000101, 2)
    1
    >>> get_bit(0b00000101, 1)
    0
    """
    return (byte & (1 << bit_num)) >> bit_num


def byte_to_bits(byte):
    """ Return the representation of a byte as a string of bits.
    @param int byte: a given byte
    @rtype: str
    >>> byte_to_bits(14)
    '00001110'
    """
    return "".join([str(get_bit(byte, bit_num))
                    for bit_num in range(7, -1, -1)])


def bits_to_byte(bits):
    """ Return int represented by bits, padded on right.
    @param str bits: a string representation of some bits
    @rtype: int
    >>> bits_to_byte("00000101")
    5
    >>> bits_to_byte("101") == 0b10100000
    True
    """
    return sum([int(bits[pos]) << (7 - pos)
                for pos in range(len(bits))])


# ====================
# Functions for compression


def make_freq_dict(text):
    """ Return a dictionary that maps each byte in text to its frequency.
    @param bytes text: a bytes object
    @rtype: dict{int,int}
    >>> d = make_freq_dict(bytes([65, 66, 67, 66]))
    >>> d == {65: 1, 66: 2, 67: 1}
    True
    """

    freq_dict = {}
    for byte in text:
        if byte in freq_dict.keys():
            freq_dict[byte] += 1
        else:
            freq_dict[byte] = 1
    return freq_dict


def huffman_tree(freq_dict):
    """ Return the root HuffmanNode of a Huffman tree corresponding
    to frequency dictionary freq_dict.
    @param dict(int,int) freq_dict: a frequency dictionary
    @rtype: HuffmanNode
    >>> freq = {2: 6, 3: 4}
    >>> t = huffman_tree(freq)
    >>> result1 = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> result2 = HuffmanNode(None, HuffmanNode(2), HuffmanNode(3))
    >>> t == result1 or t == result2
    True
    """

    if len(freq_dict) < 2:
        raise ValueError('Huffman compression not suitable.')

    # freq_dict copy - for creation of tree
    temp_freq = {}
    for item in [freq_dict.items()]:
        for key, value in item:
            temp_freq[key] = value

    # Creates a set of all leafs in tree
    leaf_symbols = set(temp_freq.keys())

    # Creates all leaf nodes
    nodes = {}
    for symbol in temp_freq:
        nodes[symbol] = HuffmanNode(symbol)

    # Set symbol for new internal node
    new_node_symb = max(leaf_symbols) + 1

    while len(nodes) > 1:
        items = sorted(temp_freq.items(), key=lambda x: x[1])

        left = items[0]
        right = items[1]

        # Create new entry with combined frequency
        temp_freq[new_node_symb] = left[1] + right[1]

        items.remove(left)
        items.remove(right)
        del temp_freq[left[0]]
        del temp_freq[right[0]]

        # Left and Right become one node
        new_node = HuffmanNode(None, nodes[left[0]], nodes[right[0]])
        del nodes[left[0]]
        del nodes[right[0]]

        nodes[new_node_symb] = new_node
        new_node_symb += 1

    # len(nodes) == 1 - nodes contains root node
    _, tree = nodes.popitem()
    return tree


def get_codes(tree):
    """ Return a dict mapping symbols from tree rooted at HuffmanNode to codes.
    @param HuffmanNode tree: a Huffman tree rooted at node 'tree'
    @rtype: dict(int,str)
    >>> tree = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> d = get_codes(tree)
    >>> d == {3: "0", 2: "1"}
    True
    """

    symb_dict = {}
    if not tree.is_leaf():
        # Postorder traversal - Updates symb_dict
        update_code_dict(tree, symb_dict, "")

    else:
        # No subsaquent nodes
        symb_dict[tree.symbol] = '0'

    return symb_dict


def number_nodes(tree):
    """ Number internal nodes in tree according to postorder traversal;
    start numbering at 0.
    @param HuffmanNode tree:  a Huffman tree rooted at node 'tree'
    @rtype: NoneType
    >>> left = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> right = HuffmanNode(None, HuffmanNode(9), HuffmanNode(10))
    >>> tree = HuffmanNode(None, left, right)
    >>> number_nodes(tree)
    >>> tree.left.number
    0
    >>> tree.right.number
    1
    >>> tree.number
    2
    """

    postorder = []
    post_order_node_get(tree, postorder)

    node_count = 0
    node_pairs = []

    for elem in postorder:
        node_pairs.append((node_count, elem))
        node_count += 1

    for num, node in node_pairs:
        node.number = num


def avg_length(tree, freq_dict):
    """ Return the number of bits per symbol required to compress text
    made of the symbols and frequencies in freq_dict, using the Huffman tree.
    @param HuffmanNode tree: a Huffman tree rooted at node 'tree'
    @param dict(int,int) freq_dict: frequency dictionary
    @rtype: float
    >>> freq = {3: 2, 2: 7, 9: 1}
    >>> left = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> right = HuffmanNode(9)
    >>> tree = HuffmanNode(None, left, right)
    >>> avg_length(tree, freq)
    1.9
    """

    new_dict = get_codes(tree)
    freq = sum([value for value in freq_dict.values()])
    chars_ = sum([freq_dict[i] * len(new_dict[i]) for i in freq_dict.keys()])
    return chars_ / freq


def generate_compressed(text, codes):
    """ Return compressed form of text, using mapping in codes for each symbol.
    @param bytes text: a bytes object
    @param dict(int,str) codes: mappings from symbols to codes
    @rtype: bytes
    >>> d = {0: "0", 1: "10", 2: "11"}
    >>> text = bytes([1, 2, 1, 0])
    >>> result = generate_compressed(text, d)
    >>> [byte_to_bits(byte) for byte in result]
    ['10111000']
    >>> text = bytes([1, 2, 1, 0, 2])
    >>> result = generate_compressed(text, d)
    >>> [byte_to_bits(byte) for byte in result]
    ['10111001', '10000000']
    """

    bit_lst = []
    texts = ''

    # Creates string rep of code
    for i in text:
        value = codes[i]
        texts = texts + value

    # Appends bits to byte list
    for i in range(0, len(texts), 8):
        bit_lst.append(bits_to_byte(texts[i:(i + 8)]))

    # Byte representation of bit_lst returned
    return bytes(bit_lst)


def tree_to_bytes(tree):
    """ Return a bytes representation of the tree rooted at tree.
    @param HuffmanNode tree: a Huffman tree rooted at node 'tree'
    @rtype: bytes
    The representation should be based on the postorder traversal of tree
    internal nodes, starting from 0.
    Precondition: tree has its nodes numbered.
    >>> tree = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> number_nodes(tree)
    >>> list(tree_to_bytes(tree))
    [0, 3, 0, 2]
    >>> left = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> right = HuffmanNode(5)
    >>> tree = HuffmanNode(None, left, right)
    >>> number_nodes(tree)
    >>> list(tree_to_bytes(tree))
    [0, 3, 0, 2, 1, 0, 0, 5]
    """

    tree_bytes = []
    post_traversal(tree, tree_bytes)
    return bytes(tree_bytes)


def num_nodes_to_bytes(tree):
    """ Return number of nodes required to represent tree (the root of a
    numbered Huffman tree).
    @param HuffmanNode tree: a Huffman tree rooted at node 'tree'
    @rtype: bytes
    """

    return bytes([tree.number + 1])


def size_to_bytes(size):
    """ Return the size as a bytes object.
    @param int size: a 32-bit integer that we want to convert to bytes
    @rtype: bytes
    >>> list(size_to_bytes(300))
    [44, 1, 0, 0]
    """

    # little-endian representation of 32-bit (4-byte)
    # int size
    return size.to_bytes(4, "little")


def compress(in_file, out_file):
    """ Compress contents of in_file and store results in out_file.
    @param str in_file: input file whose contents we want to compress
    @param str out_file: output file, where we store our compressed result
    @rtype: NoneType
    """

    with open(in_file, "rb") as f1:
        text = f1.read()
    freq = make_freq_dict(text)
    tree = huffman_tree(freq)
    codes = get_codes(tree)
    number_nodes(tree)
    print("Bits per symbol:", avg_length(tree, freq))
    result = (num_nodes_to_bytes(tree) + tree_to_bytes(tree) +
              size_to_bytes(len(text)))
    result += generate_compressed(text, codes)
    with open(out_file, "wb") as f2:
        f2.write(result)


# ====================
# Functions for decompression

def generate_tree_general(node_lst, root_index):
    """ Return the root of the Huffman tree corresponding
    to node_lst[root_index].
    The function assumes nothing about the order of the nodes in the list.
    @param list[ReadNode] node_lst: a list of ReadNode objects
    @param int root_index: index in the node list
    @rtype: HuffmanNode
    #>>> lst = [ReadNode(1, 0, 1, 0), ReadNode(0, 5, 0, 7), \
    #ReadNode(0, 10, 0, 12)]
    #>>> generate_tree_general(lst, 0)
    >>> lst = [ReadNode(0, 5, 0, 7), ReadNode(0, 10, 0, 12), \
    ReadNode(1, 1, 1, 0)]
    >>> generate_tree_general(lst, 2)
    HuffmanNode(None, HuffmanNode(None, HuffmanNode(10, None, None), \
HuffmanNode(12, None, None)), \
HuffmanNode(None, HuffmanNode(5, None, None), HuffmanNode(7, None, None)))
    """

    nodes = {}

    for i in range(len(node_lst)):
        node = HuffmanNode(None)
        node.number = i
        nodes[i] = node

    # Sequences left and right 'branches' of Huffman Tree
    for i in range(len(node_lst)):

        if node_lst[i].l_type == 0:
            nodes[i].left = HuffmanNode(node_lst[i].l_data)
        else:
            nodes[i].left = nodes[node_lst[i].l_data]

        if node_lst[i].r_type == 0:
            nodes[i].right = HuffmanNode(node_lst[i].r_data)
        else:
            nodes[i].right = nodes[node_lst[i].r_data]

    # Root of tree becomes dict value at root_index
    tree = nodes[root_index]
    return tree


def generate_tree_postorder(node_lst, root_index):
    """ Return the root of the Huffman tree corresponding
    to node_lst[root_index].
    The function assumes that the list represents a tree in postorder.
    @param list[ReadNode] node_lst: a list of ReadNode objects
    @param int root_index: index in the node list
    @rtype: HuffmanNode
    >>> lst = [ReadNode(0, 5, 0, 7), ReadNode(0, 10, 0, 12), \
ReadNode(1, 0, 1, 0)]
    >>> generate_tree_postorder(lst, 2)
    HuffmanNode(None, HuffmanNode(None, HuffmanNode(5, None, None), \
HuffmanNode(7, None, None)), \
HuffmanNode(None, HuffmanNode(10, None, None), HuffmanNode(12, None, None)))
    """

    tree = None
    nodes = []
    index = 0

    for item in node_lst:
        node = HuffmanNode(None)

        # Creates left node
        if item.l_type == 0:
            node.left = HuffmanNode(item.l_data)
        else:
            node.left = nodes.pop(0)

        # Creates right node
        if item.r_type == 0:
            node.right = HuffmanNode(item.r_data)
        else:
            node.right = nodes.pop(0)

        nodes.append(node)

        if index == root_index:
            tree = node

        index += 1
    return tree


def generate_uncompressed(tree, text, size):
    """ Returns uncompresed text usin Huffman tree, tree,
    to decompress size bytes from text.
    @param HuffmanNode tree: a HuffmanNode tree rooted at 'tree'
    @param bytes text: text to decompress
    @param int size: how many bytes to decompress from text.
    @rtype: bytes
    """

    decompress_text = []
    index_bit = 0
    index_byte = 1
    bits = byte_to_bits(text[0])  # Retrieve bits of first byte

    while size:
        current_node = tree  # Starts at front
        while not current_node.is_leaf():
            if index_bit == 8:
                bits = byte_to_bits(text[index_byte])
                index_byte += 1
                index_bit = 0

            bit = bits[index_bit]
            index_bit += 1

            if bit == '0':  # '0' for left move by Huffman Algorithm
                current_node = current_node.left
            else:
                current_node = current_node.right

        decompress_text.append(current_node.symbol)
        size -= 1

    return bytes(decompress_text)


def bytes_to_nodes(buf):
    """ Return a list of ReadNodes corresponding to the bytes in buf.
    @param bytes buf: a bytes object
    @rtype: list[ReadNode]
    >>> bytes_to_nodes(bytes([0, 1, 0, 2]))
    [ReadNode(0, 1, 0, 2)]
    """

    lst = []
    for i in range(0, len(buf), 4):
        l_type = buf[i]
        l_data = buf[i+1]
        r_type = buf[i+2]
        r_data = buf[i+3]
        lst.append(ReadNode(l_type, l_data, r_type, r_data))
    return lst


def bytes_to_size(buf):
    """ Return the size corresponding to the
    given 4-byte little-endian representation.
    @param bytes buf: a bytes object
    @rtype: int
    >>> bytes_to_size(bytes([44, 1, 0, 0]))
    300
    """

    return int.from_bytes(buf, "little")


def uncompress(in_file, out_file):
    """ Uncompress contents of in_file and store results in out_file.
    @param str in_file: input file to uncompress
    @param str out_file: output file that will hold the uncompressed results
    @rtype: NoneType
    """

    with open(in_file, "rb") as f:
        num_nodes = f.read(1)[0]
        buf = f.read(num_nodes * 4)
        node_lst = bytes_to_nodes(buf)
        tree = generate_tree_general(node_lst, num_nodes - 1)
        size = bytes_to_size(f.read(4))
        with open(out_file, "wb") as g:
            text = f.read()
            g.write(generate_uncompressed(tree, text, size))

# =========================================
# Other Functions


def improve_tree(tree, freq_dict):
    """ Improve the tree as much as possible, without changing its shape,
    by swapping nodes. The improvements are with respect to freq_dict.
    @param HuffmanNode tree: Huffman tree rooted at 'tree'
    @param dict(int,int) freq_dict: frequency dictionary
    @rtype: NoneType
    >>> left = HuffmanNode(None, HuffmanNode(99), HuffmanNode(100))
    >>> right = HuffmanNode(None, HuffmanNode(101), \
    HuffmanNode(None, HuffmanNode(97), HuffmanNode(98)))
    >>> tree = HuffmanNode(None, left, right)
    >>> freq = {97: 26, 98: 23, 99: 20, 100: 16, 101: 15}
    >>> improve_tree(tree, freq)
    >>> avg_length(tree, freq)
    2.31
    """

    nodes_list = {}
    # Generate a list of symbols sorted by their frequencies
    freq_lst = sorted(freq_dict.items(), key=lambda x: x[1])

    traverse_tree(tree, nodes_list, 0)
    # Generate Tupples - Node paired with respective depth
    # Sorted by depth in descending order
    nodes = sorted(list(nodes_list.values()), key=lambda x: x[0], reverse=True)

    # Greater depth must get lesser frequency
    for i in range(len(freq_lst)):
        nodes[i][1].symbol = freq_lst[i][0]

# =========================================
# My Helper Functions


def update_code_dict(tree, symb_dict, current_code):
    """ Updates symbol dictionary, symb_dict, with elements of tree as keys and
    string, current code, as values.
    @param tree: A Huffman Tree
    @param symb_dict:  Dictionary of values and codes
    @param current_code: Str - Subsaquent code afte traversal
    @rtype: None
    >>> symb_dict = {}
    >>> tree = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> update_code_dict(tree, symb_dict, "")
    >>> print(symb_dict)
    {3: '0', 2: '1'}
    """

    if not tree.is_leaf():
        update_code_dict(tree.left, symb_dict, current_code + "0")
        update_code_dict(tree.right, symb_dict, current_code + "1")
    else:
        symb_dict[tree.symbol] = current_code


def post_order_node_get(tree, node_list):
    """ Appends all nodes in Huffman tree to node_list, through post order
    traversal.
    @param tree: An instance of a Huffman Tree.
    @param node_list: An updated list of nodes
    @rtype: None
    >>> node_list = []
    >>> left = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> right = HuffmanNode(None, HuffmanNode(9), HuffmanNode(10))
    >>> tree = HuffmanNode(None, left, right)
    >>> post_order_node_get(tree, node_list)
    >>> print(node_list)
    [HuffmanNode(None, HuffmanNode(3, None, None), HuffmanNode(2, None, None))\
, HuffmanNode(None, HuffmanNode(9, None, None), HuffmanNode(10, None, None)), \
HuffmanNode(None, HuffmanNode(None, HuffmanNode(3, None, None), \
HuffmanNode(2, None, None)), HuffmanNode(None, HuffmanNode(9, None, None), \
HuffmanNode(10, None, None)))]
    """

    if tree.is_leaf() or tree is None:
        return
    else:
        post_order_node_get(tree.left, node_list)
        post_order_node_get(tree.right, node_list)
        node_list.append(tree)


def traverse_tree(tree, nodes_list, depth):
    """Updated dict nodes_list, with nodes as key, and their depth as value.
    @param tree: A huffman tree
    @param nodes_list: The depth of the leaf node
    @param depth: current depth of the path
    @rtype: None
    >>> tree = HuffmanNode(None, HuffmanNode(101), \
    HuffmanNode(None, HuffmanNode(97), HuffmanNode(98)))
    >>> nodes_list = {}
    >>> depth = 0
    >>> traverse_tree(tree, nodes_list, depth)
    >>> print(nodes_list)
    {101: (1, HuffmanNode(101, None, None)), \
97: (2, HuffmanNode(97, None, None)), 98: (2, HuffmanNode(98, None, None))}
    """

    if tree is None:
        return
    elif tree.is_leaf():
        nodes_list[tree.symbol] = (depth, tree)
    else:
        traverse_tree(tree.left, nodes_list, depth + 1)
        traverse_tree(tree.right, nodes_list, depth + 1)


def post_traversal(huff_tree, bytes_list):
    """Append byte representation of a Huffman tree to list tree_bytes,
    traversed in post order.
    @param huff_tree:  A Huffman Tree
    @param bytes_list: List representation of byes in tree.
    @rtype: None
    >>> tree = HuffmanNode(None, HuffmanNode(3), HuffmanNode(2))
    >>> tree_bytes = []
    >>> post_traversal(tree, tree_bytes)
    >>> print(tree_bytes)
    [0, 3, 0, 2]
    """

    if huff_tree is None:
        return
    if huff_tree.is_leaf():
        return
    else:
        post_traversal(huff_tree.left, bytes_list)
        post_traversal(huff_tree.right, bytes_list)
        if huff_tree.left.is_leaf():
            bytes_list.append(0)
            bytes_list.append(huff_tree.left.symbol)
        else:
            bytes_list.append(1)
            bytes_list.append(huff_tree.left.number)

        if huff_tree.right.is_leaf():
            bytes_list.append(0)
            bytes_list.append(huff_tree.right.symbol)
        else:
            bytes_list.append(1)
            bytes_list.append(huff_tree.right.number)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import time

    mode = input("Press c to compress or u to uncompress: ")
    if mode == "c":
        fname = input("File to compress: ")
        start = time.time()
        compress(fname, fname + ".huf")
        print("compressed {} in {} seconds."
              .format(fname, time.time() - start))
    elif mode == "u":
        fname = input("File to uncompress: ")
        start = time.time()
        uncompress(fname, fname + ".orig")
        print("uncompressed {} in {} seconds."
              .format(fname, time.time() - start))