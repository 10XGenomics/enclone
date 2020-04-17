# How to compile the `enclone` software (for experts!)

You should only have to do this if you want to experimentally modify the `enclone` software.

1. For now, you can run on an x86-64 linux server or a Mac.

2. You need to have the Rust compiler installed. Detailed instructions on how to do this
can be found [here](https://www.rust-lang.org/tools/install). You can confirm that you 
have successfully installed the Rust compiler by running `rustc --version`.

3. Clone the `enclone` repository and build `enclone` using Cargo (which comes with Rust) by running:
```
git clone git@github.com:10XGenomics/enclone.git
cargo build --release --bin enclone
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; and then put `target/release/enclone` in your path. This can be done in multiple ways, including running
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`PATH=$PATH:~/where/you/cloned/enclone/target/release/enclone`. 
See also [this helpful post on StackExchange](https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Compilation takes 8-10 minutes on a 2017 MacBook Pro with a dual-core i7 and 5-7 minutes on a similar Linux machine. 

4. Copy the directory `enclone/test/inputs` to somewhere you can point to, or just leave it 
where it is.  These are test data you can play with; you can also supply your own output
from a Cell Ranger immune profiling run (so long as there is an `all_contig_annotations.json` output). 
When you read the documentation at step 6, you'll get to a place where you put `PRE=enclone/test/inputs` 
or instead with the path where your copied data reside.  But you need to supply `PRE` with a path that 
makes sense relative to your working directory.

5. Type `enclone help`, and read the terminal setup instructions there.

6. Type `enclone` and study the documentation shown there.

7. If you want to run the built-in tests, type
```
cargo test --release -- --nocapture
```

If you have problems, please write to us at enclone@10xgenomics.com.
