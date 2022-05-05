f = open('QUERY.fasta', 'r')
lines = f.readlines()
n = int(len(lines))

print('create table query (id bigserial, seq aa_seq);')
print('INSERT INTO query (seq) VALUES')
for i in range(n):
    if i % 2 == 1:
        print('(\'' + lines[i][:-1] + '\')' + (',' if i < n - 1 else ';'))
