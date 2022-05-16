#!/bin/python
import csv
import heapq
import sys
import operator

IGNORED_COLUMNS = 1
DEFAULT_REQUIRED_SEQUENCES = 100

import sys
if sys.version_info[0] > 2:
    # Python3 ?
    xrange = range
    new_open = open
    def old_open(filename, mode):
        if 'b' in mode:
            return new_open(filename, mode.replace('b', ''), newline = '')
        else:
            return new_open(filename, mode)
    open = old_open


class PrefixCell(object):
    def __init__(self, previous_cell = None, letter = '', likelihood = float('-inf')):
        self.previous_cell = previous_cell
        self.letter = letter
        self.likelihood = likelihood

    def iter_cells(self):
        cell = self
        while cell is not None:
            yield cell
            cell = cell.previous_cell

    def full_prefix(self):
        # skips the last cell (which has previous_cell == None)
        return [cell.letter for cell in self.iter_cells()][-2::-1]

    def __lt__(self, other):
        """
        Python3.* uses __lt__
        """
        if isinstance(other, PrefixCell):
            return self.likelihood < other.likelihood
        return super(PrefixCell, self) < other
            
    def __cmp__(self, other):
        """
        Python2.* uses __cmp__
        """
        return (other < self) - (self < other)

def find_most_likely_sequences(letters, rows, required_sequences):
    """
    This is the main calculation.
    """
    prefix_cells = ([PrefixCell()] * (len(letters) - 1)) + [PrefixCell(likelihood = 0)]
    for row in rows:
        new_prefixes = [[PrefixCell(previous_cell = previous_cell,
                                    letter = letter,
                                    likelihood = previous_cell.likelihood + letter_likelihood)
                         for previous_cell in prefix_cells]
                        for letter, letter_likelihood in zip(letters, list(row)[IGNORED_COLUMNS:])]
        prefix_cells = list(heapq.merge(*new_prefixes))[-required_sequences:]
    return [(prefix_cell.full_prefix(), prefix_cell.likelihood)
            for prefix_cell in prefix_cells][::-1] # reverse order - show most likely first.


def main(file_obj, required_sequences, output_filename, output_format):
    reader = csv.reader(file_obj)
    letters = next(reader)[IGNORED_COLUMNS:]
    # assert all(len(letter) == 1 for letter in letters), "Invalid letter was found in first row."
    
    sequences_likelihoods = find_most_likely_sequences(letters,
                                                       [map(float, row) for row in reader],
                                                       required_sequences)
    out = sys.stdout
    if output_format == 'csv':
        if output_filename is not None:
            out = open(output_filename, 'wb')
        writer = csv.writer(out)
        for sequence, likelihood in sequences_likelihoods:
            writer.writerow([str(likelihood)] + list(sequence))
    elif output_format == 'txt':
        if output_filename is not None:
            out = open(output_filename, 'wb')        
        for index, (sequence, likelihood) in enumerate(sequences_likelihoods):
            out.write(">%d_%f\n" % (index + 1, likelihood))
            out.write(''.join(sequence) + '\n')
    if out is not sys.stdout:
        out.close()
        

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(description = "Finds the most likely sequences")
    parser.add_option("-i", "--file", dest = "input_filename", help = "input CSV file (default stdin)", metavar="FILE")
    parser.add_option("-o", "--output", dest = "output_filename", help = "output filename (default stdout)", metavar = "FILE")
    parser.add_option("-k", "--required", dest = "required_sequences", type="int", help = "required sequences (K) (default: %d)" % (DEFAULT_REQUIRED_SEQUENCES,), default = DEFAULT_REQUIRED_SEQUENCES)
    parser.add_option("-f", "--format", dest = "output_format", help = "output format (default: txt)", type = 'choice', choices = ("txt", "csv"), default = "txt")

    (options, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("Unexpected args")

    if options.input_filename is None:
        import warnings
        warnings.warn("Missing input filename - using stdin")
        input_file_obj = sys.stdin
    else:
        input_file_obj = open(options.input_filename,'rb')
    main(input_file_obj, options.required_sequences, options.output_filename, options.output_format)
        
        
        
    
