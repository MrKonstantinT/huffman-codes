//! Simple library which provides tools that can be used to generate or read Huffman codes.
//!
//! Some areas of concern are the runtime of decoding and memory
//! required to store the Huffman tree; the latter is beneficial if we would require to index into
//! the compressed data in RAM and there is a shortage of RAM.
//!
//! # Usage
//!
//! This crate is not published on `crates.io` and can be used by adding `huffman_codes` under the
//! `dependencies` section name in your project's `Cargo.toml` as follows:
//!
//! ```toml
//! [dependencies]
//! huffman_codes = { git = "https://github.com/MrKonstantinT/huffman-codes" }
//! ```
//!
//! and the following to your crate root:
//!
//! ```rust
//! extern crate huffman_codes;
//! ```
//!
//! # Examples
//!
//! The following example shows how we can compress a sequence of chars.
//!
//! ```{.rust}
//! extern crate huffman_codes;
//!
//! use huffman_codes::{BitArrayList, MinBinaryHeap, HuffmanTree};
//!
//! fn main() {
//!     let sequence = "ccccaaaaabbbbcccccbbb".to_string();
//!     let alpha_freq = vec![('a', 5), ('b', 7), ('c', 9)];
//!     // Initialise queue.
//!     let mut queue = MinBinaryHeap::new();
//!     for &(symbol, frequency) in alpha_freq.iter() {
//!         queue.insert(HuffmanTree::new_leaf(symbol.clone(), frequency));
//!     }
//!     // Construct Huffman codes.
//!     let huffman_tree_root = HuffmanTree::new(&mut queue);
//!     // Get Huffman codes.
//!     let huffman_codes = HuffmanTree::get_codes(&huffman_tree_root);
//!     // Encode user input.
//!     let mut encoded = BitArrayList::new();
//!     for symbol in sequence.chars() {
//!         encoded.concatenate(huffman_codes.get(&symbol).unwrap().clone());
//!     }
//!     // Compressed string takes:
//!     assert_eq!(encoded.bytes().len(), 5);
//! }
//! ```
//!
//! And this example shows how we can expand compressed data to get original sequence:
//!
//! ```{.rust}
//! extern crate huffman_codes;
//!
//! use huffman_codes::{BitArrayList, MinBinaryHeap, HuffmanTree};
//!
//! fn main() {
//!     let c_sequence = BitArrayList::from(vec![10, 171, 252, 31, 128], 33);
//!     let alpha_freq = vec![('a', 5), ('b', 7), ('c', 9)];
//!     // Initialise queue.
//!     let mut queue = MinBinaryHeap::new();
//!     for &(symbol, frequency) in alpha_freq.iter() {
//!         queue.insert(HuffmanTree::new_leaf(symbol.clone(), frequency));
//!     }
//!     // Construct Huffman codes.
//!     let huffman_tree_root = HuffmanTree::new(&mut queue);
//!     // Decode bits to original input which created them.
//!     let mut buffer = Vec::new();
//!     HuffmanTree::traverse(&huffman_tree_root, &c_sequence, &mut buffer);
//!     let sequence: String = buffer.into_iter().collect();
//!     // Compressed string takes:
//!     assert_eq!(sequence, "ccccaaaaabbbbcccccbbb");
//! }
//! ```
extern crate bit_array_list;
extern crate min_binary_heap;


use std::hash::Hash;
use std::cmp::Ordering;
use std::collections::HashMap;

pub use bit_array_list::BitArrayList;
pub use min_binary_heap::MinBinaryHeap;


/// A nested binary tree type over type of symbols to compress.
///
/// Leaf nodes contain a symbol from the alphabet to compress and the frequency of that symbol.
/// Internal nodes only have the sum of the frequency of their immediate childen as their node
/// label. The type exports three main functions. One to create the tree (returns you tree root),
/// one to get the variable length codes from the tree given when you provide a tree root and
/// finally one to decompress data given you provide a tree root again and a buffer for expanded
/// data.
#[derive(Debug, Eq)]
pub enum HuffmanTree<T>
    where T: Clone + Hash + Ord
{
    InternalNode {
        freq: usize,
        left_c: Box<HuffmanTree<T>>,
        right_c: Box<HuffmanTree<T>>,
    },
    LeafSymbol { symbol: T, freq: usize },
}

impl<T> HuffmanTree<T>
    where T: Clone + Hash + Ord
{
    /// Returns a map of symbols and variable length codes.
    ///
    /// Requires a the root of a `HuffmanTree` to work.
    ///
    /// # Examples
    ///
    /// Basic usage:
    ///
    /// ```
    /// use huffman_codes::{MinBinaryHeap, HuffmanTree};
    /// let alpha_freq = vec![('a', 5), ('b', 7), ('c', 9)];
    /// // Initialise queue.
    /// let mut queue = MinBinaryHeap::new();
    /// for &(symbol, frequency) in alpha_freq.iter() {
    ///     queue.insert(HuffmanTree::new_leaf(symbol.clone(), frequency));
    /// }
    /// // Construct Huffman codes.
    /// let huffman_tree_root = HuffmanTree::new(&mut queue);
    /// // Get Huffman codes.
    /// let huffman_codes = HuffmanTree::get_codes(&huffman_tree_root);
    /// assert_eq!(huffman_codes.get(&'a').unwrap().bytes(), &vec![128]);
    /// assert_eq!(huffman_codes.get(&'b').unwrap().bytes(), &vec![192]);
    /// assert_eq!(huffman_codes.get(&'c').unwrap().bytes(), &vec![0]);
    /// ```
    pub fn get_codes(root: &HuffmanTree<T>) -> HashMap<T, BitArrayList> {
        match root {
            &HuffmanTree::LeafSymbol { .. } => panic!("Attempted 'get_codes' with a leaf node."),
            &HuffmanTree::InternalNode { left_c: ref lc, right_c: ref rc, .. } => {
                // Set up for DFS.
                let mut huffman_codes = HashMap::new();
                let mut path = BitArrayList::new();
                // Perform DFS on left child of root first.
                path.push(0);
                dfs(lc, &mut path, &mut huffman_codes);
                // Update path for left child return.
                path.pop();
                //Then DFS on right child of root.
                path.push(1);
                dfs(rc, &mut path, &mut huffman_codes);

                huffman_codes
            }
        }
    }

    /// Helper method to read the frequency label of `HuffmanTree` node.
    ///
    /// Note that this method requires a `HuffmanTree` node to be be un-`Box`ed.
    fn freq(&self) -> usize {
        match *self {
            HuffmanTree::InternalNode { freq: f, .. } => f,
            HuffmanTree::LeafSymbol { freq: f, .. } => f,
        }
    }

    /// Decode the compressed information into a buffer.
    ///
    /// Requires a `HuffmanTree` root node to depict the alphabet variable length codes.
    ///
    /// # Remarks
    ///
    /// Better ways exist to decode compressed sequence e.g. using lookup tables. This simple
    /// method just reads the information bit by bit---each bit tells us what path to take to
    /// reach a leaf node which tells us which code we just read---and it seems to be the best
    /// option for a small memory footprint while decoding (lookup tables can grow [very large][1]
    /// just to store small amounts of codes).
    ///
    /// [1]: http://commandlinefanatic.com/cgi-bin/showarticle.cgi?article=art007
    ///
    /// # Examples
    ///
    /// Basic usage:
    ///
    /// ```
    /// use huffman_codes::{BitArrayList, MinBinaryHeap, HuffmanTree};
    /// let c_sequence = BitArrayList::from(vec![15, 255, 245, 64, 171, 96], 44);
    /// let alpha_freq = vec![('a', 5), ('b', 7), ('c', 9), ('d', 2)];
    /// // Initialise queue.
    /// let mut queue = MinBinaryHeap::new();
    /// for &(symbol, frequency) in alpha_freq.iter() {
    ///     queue.insert(HuffmanTree::new_leaf(symbol.clone(), frequency));
    /// }
    /// // Construct Huffman codes.
    /// let huffman_tree_root = HuffmanTree::new(&mut queue);
    /// // Decode bits to original input which created them.
    /// let mut buffer = Vec::new();
    /// HuffmanTree::traverse(&huffman_tree_root, &c_sequence, &mut buffer);
    /// let sequence: String = buffer.into_iter().collect();
    /// // Compressed string takes:
    /// assert_eq!(sequence, "ccccaaaaabbbbcccccbbbdd");
    /// ```
    pub fn traverse(tree_root: &HuffmanTree<T>, directions: &BitArrayList, buf: &mut Vec<T>) {
        // Traverse tree in order given by compressed information.
        let mut current_node = tree_root;
        for i in 0..directions.len() {
            // Work out which is the next node.
            let next_node = match directions.is_set(i) {
                false => current_node.child(true).unwrap(),
                true => current_node.child(false).unwrap(),
            };
            // Go to next node.
            current_node = match next_node {
                &HuffmanTree::LeafSymbol { symbol: ref s, .. } => {
                    // Read an encoded symbol. Decode it.
                    buf.push(s.clone());
                    // Start back at root.
                    tree_root
                }
                internal_node => internal_node,
            };
        }
    }

    /// Perform the greedy algorithm Huffman on the alphabet and its frequency statistics to
    /// generate a `HuffmanTree`.
    ///
    /// Call this method with a min-priority queue initialised with the leaf nodes of the Huffman
    /// tree.
    ///
    /// # Remarks
    ///
    /// An important observation is that the Huffman algorithm is correct so whatever order the
    /// queue gets elements out---the order may change for symbols from alphabet with equal
    /// frequencies depending on the order inserted---the resulting Huffman tree will produce
    /// complety valid codes that vary. The problem in that is if the order varies, e.g. iterating
    /// on a `HashMap` to initialise the queue that was used to compute the alphabet histogram, we
    /// will not be able to always decode the compressed information (we may have computed
    /// different codes!). It is imperative that if varying codes are possible the order of
    /// inserting leaf nodes in the queue is always the same. This issue is prominent if the symbol
    /// probability distribution is close to uniform.
    ///
    /// This Huffman implementation uses a min-priority queue although a better solution
    /// exists where the heap is replaced by a van Boas tree which improves the runtime to
    /// _O(nlog(log(n)))_.
    ///
    /// # Examples
    ///
    /// Basic usage:
    ///
    /// ```
    /// use huffman_codes::{MinBinaryHeap, HuffmanTree};
    /// let alpha_freq = vec![('a', 5), ('b', 7), ('c', 9)];
    /// // Initialise queue.
    /// let mut queue = MinBinaryHeap::new();
    /// for &(symbol, frequency) in alpha_freq.iter() {
    ///     queue.insert(HuffmanTree::new_leaf(symbol.clone(), frequency));
    /// }
    /// // Construct Huffman codes.
    /// let huffman_tree_root = HuffmanTree::new(&mut queue);
    /// let root_label = match huffman_tree_root {
    ///     HuffmanTree::InternalNode { freq: f, .. } => f,
    ///     HuffmanTree::LeafSymbol { .. } => panic!("Expected root to be a internal node."),
    /// };
    ///
    /// assert_eq!(root_label, 21);
    /// ```
    pub fn new(queue: &mut MinBinaryHeap<HuffmanTree<T>>) -> HuffmanTree<T> {
        // Do number of symbols - 1 merging operation to create final Huffman tree.
        for _ in 0..queue.size() - 1 {
            // Merge least frequent labelled Huffman trees in queue together. Children of node z.
            let left_c = queue.extract_min().unwrap();
            let right_c = queue.extract_min().unwrap();
            // Allocate new node z with label = the sum of the frequencies of its children.
            let node_z = new_int_node(left_c.freq() + right_c.freq(), left_c, right_c);
            // Find next least frequent labelled Huffman trees to merge.
            queue.insert(node_z);
        }
        // Return root of Huffman tree.
        queue.extract_min().unwrap()
    }

    /// Helper method that creates a `HuffmanTree` leaf node from a symbol its frequencies.
    ///
    /// # Examples
    ///
    /// Basic usage:
    ///
    /// ```
    /// use huffman_codes::HuffmanTree;
    /// let leaf = HuffmanTree::new_leaf('a', 32);
    ///
    /// assert_eq!(leaf, HuffmanTree::LeafSymbol { symbol: 'a', freq: 32 });
    /// ```
    pub fn new_leaf(s: T, f: usize) -> HuffmanTree<T> {
        HuffmanTree::LeafSymbol {
            symbol: s,
            freq: f,
        }
    }
}

impl<T> Ord for HuffmanTree<T>
    where T: Clone + Hash + Ord
{
    fn cmp(&self, other: &HuffmanTree<T>) -> Ordering {
        self.freq().cmp(&other.freq())
    }
}

impl<T> PartialOrd for HuffmanTree<T>
    where T: Clone + Hash + Ord
{
    fn partial_cmp(&self, other: &HuffmanTree<T>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> PartialEq for HuffmanTree<T>
    where T: Clone + Hash + Ord
{
    fn eq(&self, other: &HuffmanTree<T>) -> bool {
        self.freq() == other.freq()
    }
}

/// Trait for getting the child node of binary tree enums which are nested.
///
/// Intention is to allow elegant `HuffmanTree` [traverse()][2] i.e. massively improve the
/// readability of the code inside that function even though in practise there is like one
/// `HuffmanTree` which is not boxed (which was no trouble when implementing DFS). This also helped
/// to create shorter tests. This is just an exploitation of function overloading.
///
/// [2]: enum.HuffmanTree.html#method.traverse
trait GetChild<T> {
    /// This method returns the required child of a node, or `None` if node does not have that
    /// child.
    ///
    /// `true` tries to select left child and `false` tries to select right child.
    fn child(&self, left_c: bool) -> Option<&T>;
}

impl<T> GetChild<HuffmanTree<T>> for HuffmanTree<T>
    where T: Clone + Hash + Ord
{
    fn child(&self, left_c: bool) -> Option<&HuffmanTree<T>> {
        let children = match self {
            &HuffmanTree::LeafSymbol { .. } => return None,
            // Cannot copy pointer as it would violate single Rust's single ownership property.
            &HuffmanTree::InternalNode { left_c: ref lc, right_c: ref rc, .. } => (lc, rc),
        };
        // Return user preference.
        match left_c {
            true => Some(children.0),
            false => Some(children.1),
        }
    }
}

impl<T> GetChild<HuffmanTree<T>> for Box<HuffmanTree<T>>
    where T: Clone + Hash + Ord
{
    fn child(&self, left_c: bool) -> Option<&HuffmanTree<T>> {
        let children = match **self {
            HuffmanTree::LeafSymbol { .. } => return None,
            // Cannot copy pointer as it would violate single Rust's single ownership property.
            HuffmanTree::InternalNode { left_c: ref lc, right_c: ref rc, .. } => (lc, rc),
        };
        // Return user preference.
        match left_c {
            true => Some(children.0),
            false => Some(children.1),
        }
    }
}

/// Helper method that creates a `HuffmanTree` internal node.
///
/// We must specify internal node's children which aligns nice with the `huffman` procedure i.e.
/// we do not require an internal node to have only one child at some point in time.
fn new_int_node<T>(f: usize, lc: HuffmanTree<T>, rc: HuffmanTree<T>) -> HuffmanTree<T>
    where T: Clone + Hash + Ord
{
    HuffmanTree::InternalNode {
        freq: f,
        left_c: Box::new(lc),
        right_c: Box::new(rc),
    }
}

/// DFS---the meat of the [get_codes()][3] function.
///
/// [3]: enum.HuffmanTree.html#method.get_codes
fn dfs<T>(huffman_tree: &Box<HuffmanTree<T>>,
          path: &mut BitArrayList,
          huffman_codes: &mut HashMap<T, BitArrayList>)
    where T: Clone + Hash + Ord
{
    match **huffman_tree {
        HuffmanTree::LeafSymbol { symbol: ref s, .. } => {
            huffman_codes.insert(s.clone(), path.clone());
        }
        HuffmanTree::InternalNode { left_c: ref lc, right_c: ref rc, .. } => {
            // Go down left child first.
            path.push(0);
            dfs(lc, path, huffman_codes);
            // Update path for left child return.
            path.pop();
            // Then go down right child.
            path.push(1);
            dfs(rc, path, huffman_codes);
            // Update path for right child return.
            path.pop();
        }
    }
}

#[cfg(test)]
mod tests {
    use HuffmanTree;
    use GetChild;

    #[test]
    fn freq_helper() {
        let (leaf1, leaf2) = setup_leafs();

        assert_eq!(leaf1.freq(), 13);
        assert_eq!(leaf2.freq(), 5);

        let internal_node = super::new_int_node(leaf1.freq() + leaf2.freq(), leaf2, leaf1);

        assert_eq!(internal_node.freq(), 18);
    }

    #[test]
    fn new_int_node_helper() {
        let (leaf1, leaf2) = setup_leafs();

        let internal_node = super::new_int_node(leaf2.freq() + leaf1.freq(), leaf2, leaf1);

        assert_eq!(internal_node.freq(), 18);
        assert_eq!(internal_node.child(true).unwrap(),
                   &HuffmanTree::LeafSymbol {
                        symbol: 'b',
                        freq: 5,
                    });
        assert_eq!(internal_node.child(false).unwrap(),
                   &HuffmanTree::LeafSymbol {
                        symbol: 'a',
                        freq: 13,
                    });
    }

    #[test]
    fn child_function_overloading() {
        let (leaf1, leaf2) = setup_leafs();
        let boxed_leaf2 = Box::new(leaf2);

        assert_eq!(leaf1.child(true), None);
        assert_eq!(leaf1.child(false), None);
        assert_eq!(boxed_leaf2.child(true), None);
        assert_eq!(boxed_leaf2.child(false), None);

        let leaf2 = *boxed_leaf2;
        let internal_node = super::new_int_node(leaf2.freq() + leaf1.freq(), leaf2, leaf1);

        assert_eq!(internal_node.child(true),
                   Some(&HuffmanTree::LeafSymbol {
                             symbol: 'b',
                             freq: 5,
                         }));
        assert_eq!(internal_node.child(false),
                   Some(&HuffmanTree::LeafSymbol {
                             symbol: 'a',
                             freq: 13,
                         }));

        let boxed_internal_node = Box::new(internal_node);

        assert_eq!(boxed_internal_node.child(true),
                   Some(&HuffmanTree::LeafSymbol {
                             symbol: 'b',
                             freq: 5,
                         }));
        assert_eq!(boxed_internal_node.child(false),
                   Some(&HuffmanTree::LeafSymbol {
                             symbol: 'a',
                             freq: 13,
                         }));
    }

    fn setup_leafs() -> (HuffmanTree<char>, HuffmanTree<char>) {
        (HuffmanTree::new_leaf('a', 13), HuffmanTree::new_leaf('b', 5))
    }
}
