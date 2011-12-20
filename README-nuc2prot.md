% Description of nuc2prot.py script
% Paul M. Magwene
% 20 December 2011

# nuc2prot.py

This module/script allows you to translate nucleotide sequences to amino acid sequences, and to "thread" nucleotide sequences onto aligned protein sequences.

Put this script somewhere in your `PATH` or in your working directory and make it executable using `chmod +x nuc2prot.py`.


## Translating 

The default action of this script is to translate a FASTA file with DNA nucleotide sequences into a corresponding FASTA file with amino acid sequences

The script can read input either from stdin or a file, and write output to either stdout or a file. Below are some examples.

- Read from stdin, write to stdout

    $ ./nuc2prot.py < gpa2.fasta 

- Read from file, write to stdout

    $ ./nuc2prot.py -i gpa2.fasta 
    
- Read from file, write to stdout, explicitly choose translation

    $ ./nuc2prot.py -i gpa2.fasta --action=translate    
     
- Read from file, write to another file

    $ ./nuc2prot.py -i gpa2.fasta -o gpa2prot.fasta
    


## Alignment

This script doesn't deal with alignment. I recommend you use [MAFFT] (http://mafft.cbrc.jp/alignment/software/) or [T-COFFEE](http://www.tcoffee.org/Projects_home_page/t_coffee_home_page.html) to align your protein sequences. If using T-COFFEE make sure the output is FASTA format.

- Running protein alignment w/MAFFT

    $ mafft gpa2prot.fasta > gpa2prot-aln.fasta


## Threading

Threading a nucleotide sequence onto a protein alignment involves two steps: 1) confirming that the codons in the DNA sequence correspond to the respective amino acids in the amino acid sequence; and 2) adding gaps where necessary to respect the alignment.

Here's how to thread a sequence onto a protein alignment:

- Read nucleotides from file, write to stdout

    $ ./nuc2prot.py -i gpa2.fasta --action=thread --protfile=gpa2prot-aln.fasta 
    
- Read nucleotides from file, write to another file

    $ ./nuc2prot.py -i gpa2.fasta -o gpa2-threaded.fasta --action=thread --protfile=gpa2prot-aln.fasta 