#!/usr/bin/python3

import sys

f = open(sys.argv[1], 'r')
lines = f.readlines()
n = int(len(lines))

print('CREATE TABLE ' + sys.argv[2] + ' (id bigserial, seq aa_seq);')
print('INSERT INTO ' + sys.argv[2] + ' (seq) VALUES')
for i in range(n):
    if i % 2 == 1:
        print('(\'' + lines[i][:-1] + '\')' + (',' if i < n - 1 else ';'))
