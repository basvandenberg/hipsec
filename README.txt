HIPSEC 2.0                                                            2013-03-11



Hipsec predicts if gene overexpression of an extracellular protein will lead 
to successful high-level production in Aspergillus niger. The prediction is 
based on a protein's amino acid sequence. More information about the prediction 
method can be found in:

Exploring sequence characteristics related to high-level production of secreted
proteins in Aspergillus niger. B.A. van den Berg, M.J.T. Reinders, M. Hulsman,
L. Wu, H.J. Pel, J.A. Roubos, D. de Ridder (2012) PLoS ONE (accepted).

The software is developed and tested in a linux environment. The instructions
in this document only apply to those who want to run the software on a linux
operating system. The python script predictor.py should also work on Mac OS and
Windows with python installed.

A webserver that runs the software is available at: 

http://helix.ewi.tudelft.nl/hipsec.



################################################################################
# SOFTWARE REQUIREMENTS
################################################################################

Either python >= 2.7 or python >= 3.2 is required to run the predictor.py
software. The PyQT4 module is only required if you want to use the GUI. On
Debian/Ubuntu style systems this module can be installed using:

sudo apt-get install python-qt4



################################################################################
# USAGE
################################################################################



### Command line application ###

The predictor.py script can be executed from the command line. It takes two
(required) arguments: -i (--input-file) to define the path to the input fasta
file that contains the protein sequence(s) (uppercase amino acid sequences),
and -o (--output_file) to define the path to the file to which the predictions
will be written. Note that if the output file already exists, this file will be
overwritten! 

./predictor.py -i input_protein_file.fsa -o predictions.txt

For prediction based on ORF codon sequences, -c or --codon_sequence should be
provided as parameter, e.g.

./predictor.py -c -i input_orf_file.fsa -o predictions



### Graphical user interface ###

A very basic GUI is provided. You can run it by double-clicking
predictor_gui.py (and choose Run). You can open a fasta file for which the
predictions will then be shown in the text area. The predictions can be saved
as txt-file or copied to the clipboard. Be careful, the software does not
check the input. No output will be generated if you open a file that is not a
fasta file or if the fasta file contains errors. An output will be given if you
load a fasta file with for example nucleotide sequences, however these results
are useless. A fasta file with (uppercase) amino acid sequences should be used.

NOTE: the GUI can only be used with protein amino acid sequences. Prediction of
ORF codon sequences is not implemented yet.



### Python module ###

The code in predictor.py can be used as python module. The module can be
imported and used as shown in the following example code:

-example.py---------------------------------------------------------------------

from hipsec import predictor

# define path to the input fasta file
fasta_f = '/path/to/some/fasta/file.fsa'

# obtain a list of predictions
predictions = predictor.get_output(fasta_f)

for (id, prediction, annotation) in predictions:
    print '%s\t%.3f\t%s' % (id, prediction, annotation)

--------------------------------------------------------------------------------

To import the predictor, the path to the hipsec directory must be in the python
path. You can set this path by adding the following line to your .bashrc file.
You need to replace /path/to/your/ with the location where you stored the
hipsec_0.1 directory.

export PYTHONPATH=$PYTHONPATH:/path/to/your/hipsec_1.0

