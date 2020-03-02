# Huffman codes

The library provides a full binary tree data structure and procedures required to construct or decode Huffman codes from the alphabet of a sequence of symbols and the associated frequency of said symbols in the sequence. With these tools lossless data compression can be achieved.

[![Build Status](https://travis-ci.org/konstantindt/huffman-codes.svg?branch=master)](https://travis-ci.org/konstantindt/huffman_codes)

### Advantages

With Huffman codes data compression is very effective: savings between `20%` to `90%` are typical (depending on the sequence of symbols) ([Cormen et. al., 2009](#references)).

I have presented a possible way of accessing the compressed symbols without decoding the whole sequence [here](https://gist.github.com/konstantindt/1cdd564ca6b286389a8a874dc642f1a6).

If the symbols in the input sequence generally occur in a non-repeating order we will not get effective compression with [Lempel-Ziv-Welch](http://cs.indstate.edu/~ngandepalli/Abstract.pdf) (LZW) which has multiple advantages including at least `50%` size reduction.

### Disadvantages

If the probability distribution of the symbols in the input sequence is close to uniform then we require a very small alphabet size to see reduction in size; for large alphabet we may start expanding the data (length of code(s) surpasses length of original symbol unit of storage) because there are not many of the most occurring symbol symbols---which have the shortest possible code---to compensate in any way.

## Documentation

To view project documentation first build it:

```
$ cargo doc
```

and then open the following file with a web browser:

```
[crate root]/target/doc/huffman_codes/index.html
```

## Usage

Add this entry under `Cargo.toml` `dependencies` section name:

```toml
[dependencies]
huffman_codes = { git = "https://github.com/konstantindt/huffman-codes" }
```

and the following to your crate root:

```rust
extern crate huffman_codes;
```

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).

## References

Cormen, T. H., Leiserson, C. E., Rivest, R. L. and Stein, C. (2009). _Introduction to Algorithms_, chapter Greedy Algorithms. 3 Ed. Massachusetts USA, The MIT Press.
