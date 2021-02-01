# utility functions for other test classes

import gzip

import numpy as np

import pytest

def _make_info_dict(info):
    d = {}
    for pair in info.split(';'):
        k, v = pair.split('=')
        vals = np.array(v.split(','))
        try:
            vals = vals.astype(float)
        except ValueError: # not a numeric field
            pass
        d[k] = vals
    return d

def _make_format_list(fmt):
    l = []
    for val in fmt.split(':'):
        vals = np.array(val.split(','))
        try:
            vals = vals.astype(float)
        except ValueError: # not a numeric field
            pass
        l.append(vals)
    return l

# fname1 should be the output file
# fname2 should be the control file
# allow reordering of header lines
def assert_same_vcf(fname1, fname2, info_ignore = set(),
                     format_ignore = set()):
    open_fn = open
    if fname1[-3:] == '.gz':
        open_fn = gzip.open
    print(fname1, fname2)
    headers1 = set()
    headers2 = set()
    with open_fn(fname1, mode='rt') as file1, open_fn(fname2, mode='rt') as file2:
        iter1 = iter(file1)
        iter2 = iter(file2)
        failed = False
        while True:
            line1, line2 = _grab_line_for_assertion(iter1, iter2)
            if line1[0] != '#':
                raise ValueError('Output VCF header truncated abruptly')
            if line1[1] != '#' and line2[1] == '#':
                print('Output VCF header has fewer lines than the control header')
                failed = True
                headers2.add(line2)
                break
            if line1[1] == '#' and line2[1] != '#':
                print('Output VCF header has more lines than the control header')
                failed = True
                headers1.add(line1)
                break
            if line1[1] != '#' and line2[1] != '#':
                break # these are the sample lines
            headers1.add(line1)
            headers2.add(line2)

        num_commands_header1 = len([0 for line in headers1 if '##command' in line])
        num_commands_header2 = len([0 for line in headers2 if '##command' in line])
        if num_commands_header1 != num_commands_header2:
            raise ValueError(
                'Found {num_commands_header1} ##command lines in the output '
                'vcf but {num_commands_header2} in the control vcf'
            )

        # ignore command lines, they will be different
        headers1 = [line for line in headers1 if not '##command' in line]
        headers2 = [line for line in headers2 if not '##command' in line]

        # compare headers
        for line in headers1:
            if line not in headers2:
                print("Found header line '" + line + "' in the output vcf "
                      "but not in the control vcf.")
                failed = True
        for line in headers2:
            if line not in headers1:
                print("Found header line '" + line + "' in the control vcf "
                      "but not in the output vcf.")
                failed = True
        if failed:
            raise ValueError("VCF headers not identical. See output")

        # compare sample lines
        if line1 != line2:
            raise ValueError('Output vcf sample line differs from control vcf'
                  ' sample line.\nSample line in output vcf'
                  ': ' + line1 + '\nSample line in control vcf: ' + line2)

        # compare loci
        iter1 = iter(file1)
        iter2 = iter(file2)
        linenum = 1 + len(headers1) + num_commands_header1
        format_ignore_idxs = set()
        while True:
            linenum += 1
            lines = _grab_line_for_assertion(iter1, iter2)
            if lines is None:
                return
            line1, line2 = lines
            line1 = line1.split()
            line2 = line2.split()
            for idx in range(len(line1)):
                if idx <= 6 or idx == 8:
                    if idx == 8:
                        fmt = line1[idx].split(':')
                        for val in format_ignore:
                            format_ignore_idxs.add(fmt.index(val))
                    if line1[idx] == line2[idx]:
                        continue
                    if idx == 0:
                        field_name = 'CHROM'
                    if idx == 1:
                        field_name = 'POS'
                    if idx == 2:
                        field_name = 'ID'
                    if idx == 3:
                        field_name = 'REF'
                    if idx == 4:
                        field_name = 'ALT'
                    if idx == 5:
                        field_name = 'QUAL'
                    if idx == 6:
                        field_name = 'FILTER'
                    if idx == 8:
                        field_name = 'FORMAT'
                    raise ValueError('Output file differs from control file'
                                     ' at line ' + str(linenum) + ' at field '
                                     + field_name + '\nOutput line: ' + line1[idx] +
                                     '\nControl line: ' + line2[idx])
                elif idx == 7:
                    # INFO field, allow permuations and changes from int to
                    # double
                    info1 = _make_info_dict(line1[7])
                    info2 = _make_info_dict(line2[7])
                    if info1.keys() != info2.keys():
                        raise ValueError(
                            'Output file differs from control file'
                             ' at line ' + str(linenum) + ' where they '
                            'have different INFO fields' + '\nOutput: ' +
                            line1[7] + '\nControl: ' + line2[7])
                    for k in info1:
                        if k in info_ignore:
                            continue
                        if not np.all(info1[k] == info2[k]):
                            raise ValueError(
                                'Output file differs from control file'
                                ' at line ' + str(linenum) + ' at INFO '
                                'field ' + k + '\nOutput: ' + str(info1[k]) +
                                '\nControl: ' + str(info2[k]))
                elif idx > 8:
                    # FORMAT field, allow changes from '.' to '.,.' and
                    # allow changes from int to double
                    format1 = _make_format_list(line1[idx])
                    format2 = _make_format_list(line2[idx])
                    sample_num = str(idx - 8)
                    if len(format1) != len(format2):
                        raise ValueError(
                            'Output file differs from control file'
                             ' at line ' + str(linenum) + ' where they '
                             'have different numbers of fields for sample #' +
                             sample_num + '\n' + line1[idx] + '\n' + line2[idx])
                    for count, (f1, f2) in enumerate(zip(format1, format2)):
                        if count in format_ignore_idxs:
                            continue
                        if (f1.dtype.kind == 'U' and
                                np.all(f1 == '.') and
                                np.all(f2 == '.')):
                            continue
                        if np.issubdtype(f1.dtype, np.floating):
                            test = pytest.approx(f1) == f2
                        else:
                            test = np.all(f1 == f2)
                        if not test:
                            raise ValueError(
                                'Output file differs from control file'
                                ' at line ' + str(linenum) + ' at sample #'
                                + sample_num + ' at field ' + str(count + 1) +
                                '\nOutput: ' + str(f1) + '\nControl: ' + str(f2))


# fname1 should be the output file
# fname2 should be the control file
# files should not be gzipped
def assert_same_file(fname1, fname2, simple_name):
    print(fname1, fname2)
    with open(fname1) as file1, open(fname2) as file2:
        iter1 = iter(file1)
        iter2 = iter(file2)
        linenum = 0
        while True:
            linenum += 1
            lines = _grab_line_for_assertion(iter1, iter2)
            if lines is None:
                return
            
            line1, line2 = lines
            if line1 != line2:
                raise ValueError(
                    'Output ' + simple_name + ' file differs from control file'
                    ' at line ' + str(linenum) + '.\nLine in output'
                    ' file: ' + line1 + '\nLine in control file: ' + line2)

# return a pair of lines
# or raise an error if one iterator returned before the other
# iter1 should be the output file
# iter2 should be the control file
def _grab_line_for_assertion(iter1, iter2):
    file1ended = False
    file2ended = False
    try:
        line1 = next(iter1)
    except StopIteration:
        file1ended = True
    try:
        line2 = next(iter2)
    except StopIteration:
        file2ended = True
    if file1ended != file2ended:
        if file1ended:
            raise ValueError(
                'Output file has fewer lines than control file. '
            )
        else:
            raise ValueError(
                'Output file has more lines than control file. '
            )
    if file1ended and file2ended:
        return None

    return line1.strip(), line2.strip()

