f = open('DB.fasta', 'r')
lines = f.readlines()
n = int(len(lines))

print('create table target (id bigint, seq aa_seq);')

for i in range(n):
    if i % 2 == 1:
        print('INSERT INTO target VALUES (' + str(int((i + 1) / 2)) + ', \'' + lines[i][:-1] + '\');')
    # print(i, ' ', lines[i])
# print(lines[0:10])
